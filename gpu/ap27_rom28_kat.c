/*
 * ap27_rom28_kat.c -- Vulkan host driver for the AP27 CRT/ROM-batch
 * sieve, KAT verification build.
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * Builds the same 32-prime ROM table as ../cpu/ap27_crt_sieve.c, uploads
 * it once, and dispatches ap27_rom28.comp across a real GPU to scan
 * forward from a known AP27 (a KAT, known-answer-test) and confirm it's
 * found at the exact predicted position -- the standard correctness
 * check used throughout this project. Real, measured throughput is
 * printed too, but that's incidental to what this file is for: proving
 * the port is correct, not benchmarking it.
 *
 * Dispatch is chunked (CHUNK_BATCHES below) rather than issued as one
 * giant submission: a single dispatch running for more than roughly a
 * couple of seconds can trip a GPU driver's own hang-detection watchdog,
 * which silently truncates the work while the Vulkan API still reports
 * success -- a real failure mode on real hardware (confirmed on both an
 * integrated AMD GPU and an ARM/V3D embedded one), not a hypothetical one.
 * completed_count (see the shader) gives a second, independent way to
 * detect this even on a chunk that happened to produce zero hits --
 * working here, including across the multi-chunk KAT `make check` runs
 * (CHUNK_BATCHES=2500000, some KATs run 2600000 batches: two chunks).
 * It is
 * DISABLED in the production driver, ap27_rom28_search.c, whose
 * per-chunk reset-and-resubmit loop reads delta=0 regardless of reset
 * strategy (see that file's own comment at the per-chunk check site);
 * a reader relying on this comment to vouch for production safety
 * would be wrong to assume it carries over. TODO: diagnose and
 * re-enable there, or replace the per-invocation global atomic with a
 * per-workgroup subgroup-reduced add -- not done yet; production's
 * real safety net today is `dmesg`, checked by hand after the fact.
 *
 * Build:   make gpu (from the repo root) -- see the Makefile for the
 *          exact flags, rather than duplicating them here where they
 *          could drift out of sync.
 * Run:     ./ap27_rom28_kat [window] [nbatches] [kat_index 1-4]
 */

#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <vulkan/vulkan.h>
#include "../common-modinv.h"
#include "../common-cliargs.h"

#define BATCH   512u
#define WORDS   16u
#define NP      28u
#define CHUNK_BATCHES 2500000u
/* must match ap27_rom28.comp's LOCAL_SIZE_X exactly -- see that
 * file's comment for the real failure mode a drifted constant causes */
#define LOCAL_SIZE_X 64u
/* must match ap27_rom28.comp's HIT_CAP exactly */
#define HIT_CAP 8192u

/* Same CRT-merge design as ../cpu/ap27_crt_sieve.c: the 7 primes
 * dividing C are CRT-merged into 3 combined-modulus groups (see the
 * README's "Merging small primes" section for the measured throughput
 * win) -- slots 0-2 are those groups, slots 3-27 are the 25 table
 * primes 37..151, unchanged. */
#define NGROUP 3u
static const uint32_t GROUP_PRIMES[NGROUP][3] = { {5,7,11}, {13,23,1}, {17,19,1} };
static const uint32_t GROUP_LEN[NGROUP] = { 3, 2, 2 };

static const uint32_t PRIMES32[NP] = {
	385, 299, 323,
	37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
	101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151
};
static const uint32_t ROM_OFFSET32[NP] = {
	0, 6160, 10944,
	16112, 16704, 17360, 18048, 18800, 19648, 20592, 21568, 22640,
	23776, 24944, 26208, 27536, 28960, 30512, 32128, 33776, 35488,
	37232, 39040, 41072, 43168, 45360, 47584, 49968
};
#define TOTAL_WORDS 52384u

static const uint64_t C_VAL = 223092870ull;

/* "Tough luck" B rejection -- see ../cpu/ap27_crt_sieve.c's own copy of
 * this function for the full mechanism/precedent writeup. Short version:
 * B % P == 0 for any of the 25 table primes (PRIMES32[NGROUP..])
 * makes that prime's table-build D identically zero, which would make
 * that whole prime's table "always bad" and permanently zero-hits the
 * rest of the run with completely normal-looking throughput. Caught
 * here first, with a clear message, rather than reaching
 * modinv_u64(0,P) at all (whose own post-condition check would abort
 * too, just less specifically). The 4 real KATs are all already
 * known-clean against this, but the check costs nothing and protects
 * any future KAT/B added here from the same silent failure mode. */
/* B%29==0 / B%31==0 is a DIFFERENT, WORSE failure mode than the
 * table-prime case above: build_good5394()'s g29/g31 arrays encode the
 * legal starting residues as {step, 2*step mod P} where step=(B*C) mod
 * P, valid only when step!=0. When B%P==0, step=0, so this formula
 * becomes {0,0} and
 * the gate ACTIVELY GENERATES ONLY A%P==0 positions -- the one genuinely
 * composite residue (all 27 terms collapse to the same value mod P when
 * step=0, so A%P==0 makes every single term divisible by P at once) --
 * as if it were legal. Since neither 29 nor 31 is among the 32 primes
 * this design's own sieve checks (both are handled entirely by this
 * gate), every resulting "survivor" sails through undetected -- not
 * reduced coverage, but a 100% silent false-positive rate with
 * completely normal-looking throughput. */
static void reject_tough_luck_b(uint64_t B)
{
	if (B % 29 == 0) {
		fprintf(stderr, "FAIL: B=%llu is divisible by 29 -- this gate's own "
			"closed form (build_good5394) degenerates when B shares this "
			"factor: it would silently generate ONLY A%%29==0 starting "
			"positions (the one genuinely composite residue, not a legal "
			"one), and since 29 isn't checked by the 32-prime sieve itself, "
			"every resulting \"survivor\" would be a silent, 100%% false "
			"positive with completely normal-looking throughput. Choose a "
			"B not divisible by 29, 31, or 37..151.\n",
			(unsigned long long)B);
		exit(1);
	}
	if (B % 31 == 0) {
		fprintf(stderr, "FAIL: B=%llu is divisible by 31 -- this gate's own "
			"closed form (build_good5394) degenerates when B shares this "
			"factor: it would silently generate ONLY A%%31==0 starting "
			"positions (the one genuinely composite residue, not a legal "
			"one), and since 31 isn't checked by the 32-prime sieve itself, "
			"every resulting \"survivor\" would be a silent, 100%% false "
			"positive with completely normal-looking throughput. Choose a "
			"B not divisible by 29, 31, or 37..151.\n",
			(unsigned long long)B);
		exit(1);
	}
	for (uint32_t p = NGROUP; p < NP; p++) {
		uint32_t P = PRIMES32[p];
		if (B % P == 0) {
			fprintf(stderr, "FAIL: B=%llu is divisible by the table prime %u -- "
				"D=0 in that prime's table-build closed form, which would make "
				"that prime's table permanently \"always bad\", zero real hits "
				"for the whole run (\"tough luck\", pick a different B)\n",
				(unsigned long long)B, P);
			exit(1);
		}
	}
}

/* Real, publicly verifiable known AP27 examples -- same set as
 * ../cpu/ap27_crt_sieve.c. */
static const uint64_t KAT_A_LIST[4] = {
	605185576317848261ull, 277699295941594831ull,
	224584605939537911ull, 973044677408627339ull,
};
static const uint64_t KAT_B_LIST[4] = {
	155368778ull, 170826477ull, 81292139ull, 174505864ull,
};
static uint64_t B;
static uint64_t KAT_A;

static void die(const char *what, VkResult r)
{
	fprintf(stderr, "FAIL: %s (VkResult=%d)\n", what, r);
	exit(1);
}

/* Distinct exit code (77, the Autotools "skip" convention) for "no
 * usable Vulkan instance/device", as opposed to a real build/runtime
 * failure -- lets the Makefile tell "this machine has no GPU" (skip)
 * apart from "the GPU port is actually broken" (fail). */
static void die_no_gpu(const char *what)
{
	fprintf(stderr, "SKIP: %s\n", what);
	exit(77);
}

/* vkMapMemory, checked -- host-visible/host-coherent mappings used
 * throughout this file essentially never fail in practice, but an
 * unchecked one that does would surface only as a confusing crash or
 * silently-wrong read much later, not here where the real cause is. */
static void map_mem(VkDevice dev, VkDeviceMemory mem, VkDeviceSize size, void **out)
{
	VkResult r = vkMapMemory(dev, mem, 0, size, 0, out);
	if (r != VK_SUCCESS) die("vkMapMemory", r);
}

static uint32_t read_file_u32(const char *path, uint32_t **out_words)
{
	FILE *f = fopen(path, "rb");
	if (!f) { perror(path); exit(1); }
	fseek(f, 0, SEEK_END);
	long sz = ftell(f);
	fseek(f, 0, SEEK_SET);
	*out_words = malloc((size_t)sz);
	if (fread(*out_words, 1, (size_t)sz, f) != (size_t)sz) { perror("fread"); exit(1); }
	fclose(f);
	return (uint32_t)sz;
}

/* Resolves ap27_rom28.spv relative to argv0's own directory, not the
 * process's current working directory -- so this binary finds its
 * shader when invoked as e.g. ./gpu/ap27_rom28_kat from the repo root,
 * not only when launched from inside gpu/ itself -- without this,
 * running from the repo root fails with a real "No such file or
 * directory", not a hypothetical one. Falls back to the plain filename
 * (CWD-relative) only
 * when argv0 has no directory component at all, e.g. found via
 * $PATH. */
static void resolve_spv_path(const char *argv0, char *out, size_t out_sz)
{
	const char *slash = strrchr(argv0, '/');
	if (!slash) {
		snprintf(out, out_sz, "ap27_rom28.spv");
		return;
	}
	snprintf(out, out_sz, "%.*s/ap27_rom28.spv", (int)(slash - argv0), argv0);
}

static uint32_t find_memory_type(VkPhysicalDevice phys, uint32_t type_bits, VkMemoryPropertyFlags want)
{
	VkPhysicalDeviceMemoryProperties mp;
	vkGetPhysicalDeviceMemoryProperties(phys, &mp);
	for (uint32_t i = 0; i < mp.memoryTypeCount; i++)
		if ((type_bits & (1u << i)) && (mp.memoryTypes[i].propertyFlags & want) == want)
			return i;
	fprintf(stderr, "FAIL: no matching memory type\n");
	exit(1);
}

static VkBuffer make_buf(VkDevice dev, VkPhysicalDevice phys, VkDeviceSize size, VkDeviceMemory *mem_out)
{
	VkResult r;
	VkBuffer buf;
	VkBufferCreateInfo bci = { .sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO,
		.size = size, .usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, .sharingMode = VK_SHARING_MODE_EXCLUSIVE };
	if ((r = vkCreateBuffer(dev, &bci, NULL, &buf)) != VK_SUCCESS) die("vkCreateBuffer", r);
	VkMemoryRequirements mr;
	vkGetBufferMemoryRequirements(dev, buf, &mr);
	VkMemoryPropertyFlags host = VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT;
	VkMemoryAllocateInfo mai = { .sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO,
		.allocationSize = mr.size, .memoryTypeIndex = find_memory_type(phys, mr.memoryTypeBits, host) };
	if ((r = vkAllocateMemory(dev, &mai, NULL, mem_out)) != VK_SUCCESS) die("vkAllocateMemory", r);
	if ((r = vkBindBufferMemory(dev, buf, *mem_out, 0)) != VK_SUCCESS) die("vkBindBufferMemory", r);
	return buf;
}

static double now_seconds(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* The mod-5394 (= 6*29*31) legal-position cycle for a given B -- see
 * ../cpu/ap27_crt_sieve.c for the full explanation. Fixed for the life
 * of a search on this B. */
static uint32_t good5394[16];
static int gap5394[16];

static void build_good5394(uint64_t b)
{
	uint64_t bc29 = (b % 29) * C_VAL % 29;
	uint64_t bc31 = (b % 31) * C_VAL % 31;
	uint8_t g29[2] = { (uint8_t)bc29, (uint8_t)((2*bc29)%29) };
	uint8_t g31[4] = { (uint8_t)bc31, (uint8_t)((2*bc31)%31), (uint8_t)((3*bc31)%31), (uint8_t)((4*bc31)%31) };
	uint32_t combined[16];
	int n = 0;
	int mod6[2] = {1,5};
	for (int mi=0; mi<2; mi++) {
		for (int gi=0; gi<2; gi++) {
			for (int hi=0; hi<4; hi++) {
				uint32_t r;
				for (r=0; r<5394; r++)
					if (r%6==(unsigned)mod6[mi] && r%29==g29[gi] && r%31==g31[hi]) break;
				combined[n++] = r;
			}
		}
	}
	for (int i=1;i<16;i++) {
		uint32_t key=combined[i];
		int j=i-1;
		while(j>=0&&combined[j]>key){combined[j+1]=combined[j];j--;}
		combined[j+1]=key;
	}
	for (int i=0;i<16;i++) good5394[i]=combined[i];
	for (int i=0;i<16;i++) gap5394[i]=(int)((combined[(i+1)%16]-combined[i]+5394)%5394);
}

/* rel_offset[s] = offset from a batch's own starting A to the s-th of
 * the next 512 legal candidate positions -- see ../cpu/ap27_crt_sieve.c. */
static uint64_t rel_offset[BATCH];

static void build_rel_offsets(void)
{
	uint64_t cum = 0; int idx=0;
	for (int block=0; block<32; block++)
		for (int g=0; g<16; g++) { rel_offset[idx++]=cum; cum += (uint64_t)gap5394[g]; }
}

int main(int argc, char **argv)
{
	const char *usage =
		"usage: ap27_rom28_kat [WINDOW=<n>] [NBATCHES=<n>] [KAT=<1..4>] (any order)\n";
	VkResult r;
	uint64_t window = 3000000ull, nbatches = 20000ull, katidx64 = 3u;
	ap27_kv_spec specs[] = {
		{ .key = "WINDOW",   .type = AP27_KV_U64, .u64_out = &window,   .min = 1, .max = 10000000000000000000ull },
		{ .key = "NBATCHES", .type = AP27_KV_U64, .u64_out = &nbatches, .min = 1, .max = UINT64_MAX },
		{ .key = "KAT",      .type = AP27_KV_U64, .u64_out = &katidx64, .min = 1, .max = 4 },
	};
	if (!ap27_parse_kv_args(argc, argv, 1, specs, 3, usage)) return 1;
	ap27_echo_kv_args("ap27_rom28_kat", specs, 3);
	uint32_t katidx = (uint32_t)katidx64;
	/* Same per-chunk dispatch/readback shape as ap27_rom28_search.c (this
	 * directory): the shader's hit_batch[] is always chunk-relative (see
	 * ap27_rom28.comp), so this driver reads hits back and resets the
	 * GPU-side hit_count after every chunk, reconstructing each hit's
	 * real absolute batch as off + hit_batch[i] with a host-side,
	 * never-wrapping uint64_t off -- the same convention the production
	 * driver uses, not a separate single-shot design. */
	KAT_A = KAT_A_LIST[katidx-1];
	B = KAT_B_LIST[katidx-1];

	reject_tough_luck_b(B);
	build_good5394(B);
	build_rel_offsets();

	/* WINDOW has no meaningful static ceiling of its own -- it has to
	 * stay well below whichever KAT's target this run picked, which
	 * isn't known until katidx is. Checked here, by name, rather than
	 * letting `KAT_A - window` underflow into a huge k below and surface
	 * only as a confusing "NOT FOUND" much later. */
	if (window + (uint64_t)good5394[0] >= KAT_A) {
		fprintf(stderr, "FAIL: WINDOW=%llu is too large for KAT #%u's target "
			"A=%llu (must stay well below it)\n",
			(unsigned long long)window, katidx, (unsigned long long)KAT_A);
		return 1;
	}

	/* Construct A0 directly as a legal position (good5394[0] plus whole
	 * 5394-cycles) rather than starting from an arbitrary integer and
	 * correcting afterward -- see ../cpu/ap27_crt_sieve.c for why this
	 * is required, not just tidier. */
	uint64_t k = (KAT_A - window - (uint64_t)good5394[0]) / 5394ull;
	uint64_t A0 = (uint64_t)good5394[0] + k * 5394ull;
	printf("KAT #%u: B=%llu, target A=%llu\n", katidx, (unsigned long long)B, (unsigned long long)KAT_A);
	printf("A0=%llu (mod 5394 = %llu)\n", (unsigned long long)A0, (unsigned long long)(A0 % 5394));

	uint32_t baseline0[NP], advanceModP[NP];
	uint32_t *romTable = calloc(TOTAL_WORDS, sizeof(uint32_t));
	uint64_t advance = 32ull * 5394ull;

	for (uint32_t p = 0; p < NP; p++) {
		uint64_t P = PRIMES32[p];
		baseline0[p] = (uint32_t)(A0 % P);
		advanceModP[p] = (uint32_t)(advance % P);
		uint32_t *tab = &romTable[ROM_OFFSET32[p]];
		int is_group = (p < NGROUP); /* the 3 CRT-merged small-prime groups */
		uint64_t Dinv = 0;
		if (!is_group) {
			uint64_t c = C_VAL % P;
			uint64_t r1 = B % P;
			uint64_t D = (P - (r1*c)%P) % P;
			Dinv = modinv_u64(D, P);
		}
		for (uint64_t base=0; base<P; base++) {
			uint32_t *mask = &tab[base*WORDS];
			for (uint32_t s=0; s<BATCH; s++) {
				uint64_t raw = (base + rel_offset[s]) % P;
				int good;
				if (is_group) {
					/* good iff raw is divisible by NONE of this group's
					 * primes -- same predicate as ../cpu/ap27_crt_sieve.c's
					 * build_rom(), just re-derived here for the GPU port. */
					good = 1;
					for (uint32_t i = 0; i < GROUP_LEN[p]; i++)
						if (raw % GROUP_PRIMES[p][i] == 0) good = 0;
				} else {
					good = ((raw * Dinv) % P) > 26;
				}
				if (good) mask[s/32] |= (1u << (s%32));
			}
		}
	}
	printf("ROM size (28 slots, merged, fixed): %u bytes (%.1f KB)\n",
		TOTAL_WORDS*4, TOTAL_WORDS*4/1024.0);

	VkApplicationInfo app = { .sType = VK_STRUCTURE_TYPE_APPLICATION_INFO,
		.pApplicationName = "ap27-crt-sieve-gpu", .apiVersion = VK_API_VERSION_1_2 };
	VkInstanceCreateInfo ici = { .sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO, .pApplicationInfo = &app };
	VkInstance instance;
	r = vkCreateInstance(&ici, NULL, &instance);
	if (r == VK_ERROR_INCOMPATIBLE_DRIVER)
		die_no_gpu("vkCreateInstance: no compatible Vulkan driver found");
	else if (r != VK_SUCCESS)
		die("vkCreateInstance", r);

	uint32_t ndev = 0;
	vkEnumeratePhysicalDevices(instance, &ndev, NULL);
	VkPhysicalDevice *devs = malloc(ndev * sizeof(*devs));
	vkEnumeratePhysicalDevices(instance, &ndev, devs);
	VkPhysicalDevice phys = VK_NULL_HANDLE;
	VkPhysicalDeviceProperties chosen = {0};
	/* AP27_GPU_INDEX overrides the automatic pick entirely -- the only
	 * way to force a specific device on a system with several, since
	 * Vulkan has no notion of a "default" GPU for us to ask for. */
	const char *idx_env = getenv("AP27_GPU_INDEX");
	if (idx_env) {
		uint32_t idx = (uint32_t)strtoul(idx_env, NULL, 10);
		if (idx >= ndev) {
			fprintf(stderr, "FAIL: AP27_GPU_INDEX=%u but only %u device(s) enumerated\n",
				idx, ndev);
			exit(1);
		}
		phys = devs[idx];
		vkGetPhysicalDeviceProperties(phys, &chosen);
	} else {
		/* Prefer a discrete GPU over an integrated one; first-enumerated
		 * wins any tie within the same class, so the automatic pick is
		 * deterministic rather than "whichever happened to enumerate
		 * last", which is what this loop used to do. */
		for (uint32_t i = 0; i < ndev; i++) {
			VkPhysicalDeviceProperties props;
			vkGetPhysicalDeviceProperties(devs[i], &props);
			int is_gpu = (props.deviceType == VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU ||
			              props.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU);
			if (!is_gpu) continue;
			int better = (phys == VK_NULL_HANDLE) ||
				(props.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU &&
				 chosen.deviceType != VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU);
			if (better) { phys = devs[i]; chosen = props; }
		}
	}
	if (phys == VK_NULL_HANDLE) die_no_gpu("no real GPU found among enumerated Vulkan devices");
	fprintf(stderr, "device: %s\n", chosen.deviceName);

	uint32_t nq = 0;
	vkGetPhysicalDeviceQueueFamilyProperties(phys, &nq, NULL);
	VkQueueFamilyProperties *qprops = malloc(nq * sizeof(*qprops));
	vkGetPhysicalDeviceQueueFamilyProperties(phys, &nq, qprops);
	uint32_t qfam = UINT32_MAX;
	for (uint32_t i = 0; i < nq; i++)
		if (qprops[i].queueFlags & VK_QUEUE_COMPUTE_BIT) { qfam = i; break; }
	if (qfam == UINT32_MAX) { fprintf(stderr, "FAIL: no compute queue family\n"); exit(1); }

	float qprio = 1.0f;
	VkDeviceQueueCreateInfo qci = { .sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO,
		.queueFamilyIndex = qfam, .queueCount = 1, .pQueuePriorities = &qprio };
	VkDeviceCreateInfo dci = { .sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO,
		.queueCreateInfoCount = 1, .pQueueCreateInfos = &qci };
	VkDevice dev;
	if ((r = vkCreateDevice(phys, &dci, NULL, &dev)) != VK_SUCCESS) die("vkCreateDevice", r);
	VkQueue queue;
	vkGetDeviceQueue(dev, qfam, 0, &queue);

	VkDeviceMemory mem_rom, mem_base, mem_adv, mem_cnt, mem_hb, mem_hp, mem_comp;
	VkBuffer buf_rom = make_buf(dev, phys, TOTAL_WORDS * sizeof(uint32_t), &mem_rom);
	VkBuffer buf_base = make_buf(dev, phys, NP * sizeof(uint32_t), &mem_base);
	VkBuffer buf_adv = make_buf(dev, phys, NP * sizeof(uint32_t), &mem_adv);
	VkBuffer buf_cnt = make_buf(dev, phys, sizeof(uint32_t), &mem_cnt);
	VkBuffer buf_hb = make_buf(dev, phys, HIT_CAP * sizeof(uint32_t), &mem_hb);
	VkBuffer buf_hp = make_buf(dev, phys, HIT_CAP * sizeof(uint32_t), &mem_hp);
	VkBuffer buf_comp = make_buf(dev, phys, sizeof(uint32_t), &mem_comp);

	void *mapped;
#define UPLOAD(mem, ptr, sz) do { map_mem(dev, mem, (sz), &mapped); memcpy(mapped, (ptr), (sz)); vkUnmapMemory(dev, mem); } while (0)
	UPLOAD(mem_rom, romTable, TOTAL_WORDS * sizeof(uint32_t));
	UPLOAD(mem_base, baseline0, NP * sizeof(uint32_t));
	UPLOAD(mem_adv, advanceModP, NP * sizeof(uint32_t));
	uint32_t zero = 0;
	UPLOAD(mem_cnt, &zero, sizeof(uint32_t));
	UPLOAD(mem_comp, &zero, sizeof(uint32_t));
#undef UPLOAD

	char spv_path[512];
	resolve_spv_path(argv[0], spv_path, sizeof(spv_path));
	uint32_t *spirv;
	uint32_t spirv_bytes = read_file_u32(spv_path, &spirv);
	VkShaderModuleCreateInfo smci = { .sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO,
		.codeSize = spirv_bytes, .pCode = spirv };
	VkShaderModule shader;
	if ((r = vkCreateShaderModule(dev, &smci, NULL, &shader)) != VK_SUCCESS) die("vkCreateShaderModule", r);

	VkDescriptorSetLayoutBinding bindings[7];
	for (int i = 0; i < 7; i++)
		bindings[i] = (VkDescriptorSetLayoutBinding){ .binding = (uint32_t)i,
			.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, .descriptorCount = 1,
			.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT };
	VkDescriptorSetLayoutCreateInfo dslci = { .sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO,
		.bindingCount = 7, .pBindings = bindings };
	VkDescriptorSetLayout dsl;
	if ((r = vkCreateDescriptorSetLayout(dev, &dslci, NULL, &dsl)) != VK_SUCCESS)
		die("vkCreateDescriptorSetLayout", r);

	VkPushConstantRange pcrange = { .stageFlags = VK_SHADER_STAGE_COMPUTE_BIT, .offset = 0, .size = sizeof(uint32_t) };
	VkPipelineLayoutCreateInfo plci = { .sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO,
		.setLayoutCount = 1, .pSetLayouts = &dsl, .pushConstantRangeCount = 1, .pPushConstantRanges = &pcrange };
	VkPipelineLayout playout;
	if ((r = vkCreatePipelineLayout(dev, &plci, NULL, &playout)) != VK_SUCCESS)
		die("vkCreatePipelineLayout", r);

	VkComputePipelineCreateInfo cpci = { .sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO,
		.stage = { .sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
		           .stage = VK_SHADER_STAGE_COMPUTE_BIT, .module = shader, .pName = "main" },
		.layout = playout };
	VkPipeline pipeline;
	if ((r = vkCreateComputePipelines(dev, VK_NULL_HANDLE, 1, &cpci, NULL, &pipeline)) != VK_SUCCESS)
		die("vkCreateComputePipelines", r);

	VkDescriptorPoolSize dps = { .type = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, .descriptorCount = 7 };
	VkDescriptorPoolCreateInfo dpci = { .sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO,
		.maxSets = 1, .poolSizeCount = 1, .pPoolSizes = &dps };
	VkDescriptorPool dpool;
	if ((r = vkCreateDescriptorPool(dev, &dpci, NULL, &dpool)) != VK_SUCCESS)
		die("vkCreateDescriptorPool", r);

	VkDescriptorSetAllocateInfo dsai = { .sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO,
		.descriptorPool = dpool, .descriptorSetCount = 1, .pSetLayouts = &dsl };
	VkDescriptorSet dset;
	if ((r = vkAllocateDescriptorSets(dev, &dsai, &dset)) != VK_SUCCESS)
		die("vkAllocateDescriptorSets", r);

	VkDescriptorBufferInfo dbi[7] = {
		{ .buffer = buf_rom, .offset = 0, .range = TOTAL_WORDS * sizeof(uint32_t) },
		{ .buffer = buf_base, .offset = 0, .range = NP * sizeof(uint32_t) },
		{ .buffer = buf_adv, .offset = 0, .range = NP * sizeof(uint32_t) },
		{ .buffer = buf_cnt, .offset = 0, .range = sizeof(uint32_t) },
		{ .buffer = buf_hb, .offset = 0, .range = HIT_CAP * sizeof(uint32_t) },
		{ .buffer = buf_hp, .offset = 0, .range = HIT_CAP * sizeof(uint32_t) },
		{ .buffer = buf_comp, .offset = 0, .range = sizeof(uint32_t) },
	};
	VkWriteDescriptorSet writes[7];
	for (int i = 0; i < 7; i++)
		writes[i] = (VkWriteDescriptorSet){ .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, .dstSet = dset,
			.dstBinding = (uint32_t)i, .descriptorCount = 1, .descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
			.pBufferInfo = &dbi[i] };
	vkUpdateDescriptorSets(dev, 7, writes, 0, NULL);

	VkCommandPoolCreateInfo cpi = { .sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO, .queueFamilyIndex = qfam };
	VkCommandPool cmdpool;
	if ((r = vkCreateCommandPool(dev, &cpi, NULL, &cmdpool)) != VK_SUCCESS)
		die("vkCreateCommandPool", r);

	VkCommandBufferAllocateInfo cbai = { .sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO,
		.commandPool = cmdpool, .level = VK_COMMAND_BUFFER_LEVEL_PRIMARY, .commandBufferCount = 1 };
	VkCommandBuffer cmd;
	if ((r = vkAllocateCommandBuffers(dev, &cbai, &cmd)) != VK_SUCCESS)
		die("vkAllocateCommandBuffers", r);

	double elapsed = 0.0;
	uint64_t total_dispatched = 0;  /* actual invocation count, >= nbatches (per-chunk ceiling rounding) */
	uint64_t total_hits = 0;
	int found = 0;
	/* Tracks the current chunk's starting A-delta from A0 -- deliberately
	 * never range-checked on every `+= chunk*advance` below (see that
	 * line); see ../cpu/ap27_crt_sieve.c's own `delta`/`delta_wrapped`
	 * comment for the full reasoning, identical here. The only place
	 * this (or delta_wrapped) is ever examined is at the per-hit
	 * reconstruction below, via an exact carry check -- not a per-chunk
	 * bound here. A0 here is a fixed, small KAT target, far below
	 * UINT64_MAX, so this is all defense-in-depth: WINDOW is already
	 * checked (above) to stay well under KAT_A. */
	uint64_t delta = 0;
	int delta_wrapped = 0;
	for (uint64_t off = 0; off < nbatches; off += CHUNK_BATCHES) {
		uint32_t chunk = (uint32_t)((nbatches - off < CHUNK_BATCHES) ? (nbatches - off) : CHUNK_BATCHES);
		uint32_t groups = (chunk + LOCAL_SIZE_X - 1) / LOCAL_SIZE_X;
		total_dispatched += (uint64_t)groups * LOCAL_SIZE_X;

		/* Fold this chunk's own batch offset into baseline0 host-side
		 * (32 cheap CPU divisions), so the shader never needs a live
		 * "batchIdx % P" -- see ap27_rom28.comp's own comment. */
		uint32_t chunkBaseline0[NP];
		for (uint32_t p = 0; p < NP; p++) {
			uint64_t P = PRIMES32[p];
			uint64_t off_term = ((uint64_t)off % P) * advanceModP[p] % P;
			chunkBaseline0[p] = (uint32_t)(((uint64_t)baseline0[p] + off_term) % P);
		}
		{
			void *mapped_cb;
			map_mem(dev, mem_base, NP * sizeof(uint32_t), &mapped_cb);
			memcpy(mapped_cb, chunkBaseline0, NP * sizeof(uint32_t));
			vkUnmapMemory(dev, mem_base);
		}

		uint32_t chunk_len = chunk;
		vkResetCommandBuffer(cmd, 0);
		VkCommandBufferBeginInfo cbbi = { .sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO };
		vkBeginCommandBuffer(cmd, &cbbi);
		/* hit_count is reset GPU-side every chunk: hit_batch[]/hit_pos[]
		 * are chunk-relative (see ap27_rom28.comp), so accumulating
		 * across chunks without resetting would let one chunk's index 5
		 * collide with another's. completed_count is NOT reset here --
		 * it's a free-running total, diffed against total_dispatched
		 * once at the very end (see below). */
		vkCmdFillBuffer(cmd, buf_cnt, 0, sizeof(uint32_t), 0);
		VkMemoryBarrier fill_barrier = { .sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER,
			.srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT,
			.dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT };
		vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
			0, 1, &fill_barrier, 0, NULL, 0, NULL);
		vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline);
		vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, playout, 0, 1, &dset, 0, NULL);
		vkCmdPushConstants(cmd, playout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(uint32_t), &chunk_len);
		vkCmdDispatch(cmd, groups, 1, 1);
		VkMemoryBarrier readback_barrier = { .sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER,
			.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT, .dstAccessMask = VK_ACCESS_HOST_READ_BIT };
		vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_HOST_BIT,
			0, 1, &readback_barrier, 0, NULL, 0, NULL);
		vkEndCommandBuffer(cmd);

		VkSubmitInfo si = { .sType = VK_STRUCTURE_TYPE_SUBMIT_INFO, .commandBufferCount = 1, .pCommandBuffers = &cmd };

		double t0 = now_seconds();
		if ((r = vkQueueSubmit(queue, 1, &si, VK_NULL_HANDLE)) != VK_SUCCESS) die("vkQueueSubmit", r);
		if ((r = vkQueueWaitIdle(queue)) != VK_SUCCESS)
			fprintf(stderr, "WARNING: vkQueueWaitIdle returned %d at offset %llu -- chunk may be truncated\n",
				r, (unsigned long long)off);
		elapsed += now_seconds() - t0;

		/* read back and reconstruct THIS chunk's hits -- HIT_CAP is
		 * sized for one chunk's worth (see ap27_rom28_search.c's own
		 * comment: ~4,459 expected per 2.5M-batch chunk), not the whole
		 * run, which is exactly why this can no longer be a single
		 * end-of-run readback once chunking is involved. */
		map_mem(dev, mem_cnt, sizeof(uint32_t), &mapped);
		uint32_t hits = *(uint32_t *)mapped;
		vkUnmapMemory(dev, mem_cnt);
		total_hits += hits;
		if (hits > HIT_CAP) {
			fprintf(stderr, "WARNING: chunk at off=%llu: hits=%u > HIT_CAP=%u -- "
				"%u hits from this chunk were dropped (shader-side bounds check)\n",
				(unsigned long long)off, hits, HIT_CAP, hits - HIT_CAP);
		}
		uint32_t nread = (hits < HIT_CAP) ? hits : HIT_CAP;
		/* hits arrive in the shader's atomic append order, which is
		 * arbitrary (not ascending abs_batch/pos) -- so a hit past the
		 * 64-bit boundary does NOT imply every later hit in this chunk
		 * is also past it, unlike the CPU driver's ascending-offset
		 * scan. Skip just the out-of-range hit and keep draining the
		 * rest of the chunk's buffer, then stop after the chunk is
		 * fully drained (checked once this chunk's hits are all read,
		 * below) rather than truncating mid-chunk and silently losing
		 * in-range hits. */
		int stop_after_chunk = 0;
		/* This chunk's reportable-offset ceiling, computed once per
		 * chunk instead of twice per hit -- see
		 * ap27_rom28_search.c's identical comment for why `room`
		 * alone can't express "even within_chunk=0 is invalid" and
		 * `chunk_unsafe` has to stay a separate check. */
		uint64_t base = 0, room = 0;
		int chunk_unsafe = delta_wrapped || __builtin_add_overflow(A0, delta, &base);
		if (!chunk_unsafe) room = UINT64_MAX - base;
		if (nread > 0) {
			map_mem(dev, mem_hb, nread * sizeof(uint32_t), &mapped);
			uint32_t hit_batch[HIT_CAP];
			memcpy(hit_batch, mapped, nread * sizeof(uint32_t));
			vkUnmapMemory(dev, mem_hb);
			map_mem(dev, mem_hp, nread * sizeof(uint32_t), &mapped);
			uint32_t hit_pos[HIT_CAP];
			memcpy(hit_pos, mapped, nread * sizeof(uint32_t));
			vkUnmapMemory(dev, mem_hp);

			for (uint32_t i = 0; i < nread; i++) {
				/* hit_batch[i] is the shader's own small, chunk-relative
				 * invocation index -- added here to the host's own
				 * wide, never-wrapping `off` to get the true absolute
				 * batch index (see ap27_rom28_search.c's identical
				 * reconstruction for why both operands must be
				 * uint64_t). */
				uint64_t abs_batch = off + (uint64_t)hit_batch[i];
				/* hit_batch[i]*advance is always safe un-checked: bounded
				 * by one chunk (CHUNK_BATCHES*advance ~4.3e11), never
				 * accumulated across chunks the way `delta` is. The only
				 * overflow check in this file, run only for an actual
				 * candidate. */
				uint64_t within_chunk = (uint64_t)hit_batch[i] * advance + rel_offset[hit_pos[i]];
				if (chunk_unsafe || within_chunk > room) {
					stop_after_chunk = 1;
					continue;
				}
				uint64_t cand_A = base + within_chunk;
				int match = (cand_A == KAT_A);
				if (match) found = 1;
				printf("SURVIVOR: A=%llu (batch=%llu, pos=%u)%s\n",
					(unsigned long long)cand_A, (unsigned long long)abs_batch, hit_pos[i],
					match ? "  <-- MATCHES KAT TARGET" : "");
			}
		}
		/* delta tracks the same A-units quantity as (off+chunk)*advance,
		 * kept incrementally (see its declaration comment) -- the
		 * per-chunk product chunk*advance is always safe un-checked
		 * (bounded, ~4.3e11), unlike its accumulation into delta across
		 * many chunks, which is exactly what the carry here catches. */
		delta_wrapped |= __builtin_add_overflow(delta, (uint64_t)chunk * advance, &delta);
		if (stop_after_chunk) {
			fprintf(stderr, "ap27_rom28_kat: stopping, A would exceed "
				"the 64-bit range (off=%llu) -- WINDOW/NBATCHES is "
				"too large for this B\n", (unsigned long long)off);
			break;
		}
	}

	map_mem(dev, mem_comp, sizeof(uint32_t), &mapped);
	uint32_t completed = *(uint32_t *)mapped;
	vkUnmapMemory(dev, mem_comp);
	if ((uint64_t)completed != total_dispatched)
		fprintf(stderr, "WARNING: completed_count=%u != total_dispatched=%llu -- dispatch "
			"truncated (a real GPU reset that still returned VK_SUCCESS), results incomplete\n",
			completed, (unsigned long long)total_dispatched);

	uint64_t candidates = (uint64_t)nbatches * BATCH;
	printf("device: %s\n", chosen.deviceName);
	printf("dispatched %llu batches (%llu candidates), elapsed=%.4fs (%.1f candidates/sec)\n",
		(unsigned long long)nbatches, (unsigned long long)candidates, elapsed, candidates / elapsed);
	printf("hits=%llu\n", (unsigned long long)total_hits);
	printf("target %s\n", found ? "FOUND" : "NOT FOUND (check window/nbatches covers the real KAT position, and that KAT_A/KAT_B match)");

	return !found;
}
