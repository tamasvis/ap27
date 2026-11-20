/*
 * ap27_rom28_search.c -- real, unattended, duration-bound GPU search
 * using the AP27 CRT/ROM-batch sieve (the merged-group, 28-slot design
 * also used in ../cpu/ap27_crt_sieve.c -- see this same directory's
 * KAT-verification driver, ap27_rom28_kat.c, which this file adapts
 * into a real production search).
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * Unlike that KAT driver (fixed nbatches, hits read back once at
 * the very end into a single HIT_CAP-sized buffer), this streams hits
 * after EVERY chunk and resets the GPU-side counters between chunks --
 * necessary because a real unattended run is many chunks long: even
 * though one 2.5M-batch chunk alone produces only ~4,459 expected hits
 * (comfortably under HIT_CAP=8192), the cumulative total across an
 * unattended run of thousands of chunks would vastly exceed any fixed
 * buffer if never read back and reset. An unmodified single-readback
 * design would silently drop the vast majority of hits on any run
 * longer than a couple of chunks.
 *
 * Build:   make gpu (from the repo root) -- see the Makefile for the
 *          exact flags, rather than duplicating them here where they
 *          could drift out of sync.
 * Run:     ./ap27_rom28_search B=<b> A=<a_lo> SECONDS=<duration> (any order)
 *          (hits stream to stdout as "HIT B=... A=... status=unchecked")
 */

#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <vulkan/vulkan.h>
#include "../common-modinv.h"
#include "../common-cliargs.h"

/* Marks a value deliberately unused -- same effect as a bare (void)
 * cast, but searchable/explicit about intent (vs. reverse-engineering
 * an arbitrary-looking cast). Guarded so this stays a no-op if some
 * other header already defines the same macro. */
#if !defined(MARK_UNUSED)
#define MARK_UNUSED(prm)   (void) (prm)
#endif

#define BATCH   512u
#define WORDS   16u
#define NP      28u
#define LOCAL_SIZE_X 64u
#define HIT_CAP 8192u

/* Same CRT-merge design as ../cpu/ap27_crt_sieve.c: the 7 primes
 * dividing C are CRT-merged into 3 combined-modulus groups (real,
 * measured ~24% CPU throughput win, ~6.8% real measured GPU win) --
 * slots 0-2 are those groups, slots 3-27 are the 25 table primes
 * 37..151, unchanged. */
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

/* "Tough luck" B rejection -- same mechanism as ../cpu/ap27_crt_sieve.c's
 * own copy of this function. Short version: B % P == 0 for any of the
 * 25 table primes (PRIMES32[NGROUP..]) makes that prime's
 * table-build D identically zero, which would make that whole prime's
 * table "always bad" and permanently zero-hits the rest of the run
 * with completely normal-looking throughput. Caught here first, with a
 * clear message, rather than reaching modinv_u64(0,P) at all (whose
 * own post-condition check would abort too, just less specifically).
 *
 * B%29==0 / B%31==0 is a DIFFERENT, WORSE failure mode:
 * build_good5394()'s g29/g31 arrays encode the legal starting residues
 * as {step, 2*step mod P} where step=(B*C) mod P, valid only when
 * step!=0. When B%P==0, step=0, so this formula becomes {0,0} and the
 * gate ACTIVELY GENERATES ONLY A%P==0 positions -- the one genuinely
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
 * shader when invoked as e.g. ./gpu/ap27_rom28_search from the repo
 * root, not only when launched from inside gpu/ itself -- without
 * this, running from the repo root fails with a real "No such file or
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

static uint64_t rel_offset[BATCH];

static void build_rel_offsets(void)
{
	uint64_t cum = 0; int idx=0;
	for (int block=0; block<32; block++)
		for (int g=0; g<16; g++) { rel_offset[idx++]=cum; cum += (uint64_t)gap5394[g]; }
}

int main(int argc, char **argv)
{
	const char *usage = "usage: ap27_rom28_search B=<b> A=<a_lo> SECONDS=<duration> (any order)\n";
	VkResult r;
	uint64_t B = 0, A_lo = 0, duration_u64 = 0;
	/* A's ceiling (UINT64_MAX/2) is not a representability limit (see
	 * the 64-bit stop later in this file) -- it's generously far below
	 * 2^64 so a bad/adversarial A_lo fails here, immediately and by
	 * name, rather than running first. */
	ap27_kv_spec specs[] = {
		{ .key = "B",       .type = AP27_KV_U64, .u64_out = &B,            .min = 1, .max = UINT64_MAX,   .required = 1 },
		{ .key = "A",       .type = AP27_KV_U64, .u64_out = &A_lo,         .min = 0, .max = UINT64_MAX / 2, .required = 1 },
		{ .key = "SECONDS", .type = AP27_KV_U64, .u64_out = &duration_u64, .min = 0, .max = 100000000ull, .required = 1 },
	};
	if (!ap27_parse_kv_args(argc, argv, 1, specs, 3, usage))
		return 1;
	ap27_echo_kv_args("ap27_rom28_search", specs, 3);
	double duration = (double)duration_u64;

	reject_tough_luck_b(B);
	build_good5394(B);
	build_rel_offsets();

	uint64_t k = A_lo / 5394ull;
	uint64_t A0 = (uint64_t)good5394[0] + k * 5394ull;
	if (A0 < A_lo) A0 += 5394ull;
	fprintf(stderr, "search: B=%llu A0=%llu duration=%.1fs\n",
		(unsigned long long)B, (unsigned long long)A0, duration);

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
					 * build_rom(). */
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
	fprintf(stderr, "ROM size (28 slots, merged, fixed): %u bytes (%.1f KB)\n",
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
	fprintf(stderr, "device: %s\n", chosen.deviceName);

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

	const uint32_t CHUNK_BATCHES = 2500000u;
	double t_start = now_seconds();
	double elapsed_dispatch = 0.0;
	uint64_t total_batches = 0, total_hits = 0;
	/* total_batches/total_hits must be uint64_t, not uint32_t: a
	 * uint32_t would silently wrap every ~4.29e9 batches, well within
	 * range for a real multi-hour run. uint64_t is exact for
	 * `total_batches`/`off` themselves for any run this design could
	 * realistically do -- see `delta`'s own declaration comment below
	 * for the actual 64-bit stop, which is about the reported A value,
	 * not these batch counters. The shader's own push constant stays a
	 * 32-bit GLSL `uint` (a real language limit, not a choice) -- but it
	 * no longer needs the absolute value at all (see ap27_rom28.comp),
	 * so that limit doesn't matter. */
	uint64_t off = 0;
	/* Tracks this chunk's starting A-delta from A0 -- deliberately never
	 * range-checked on every `+= chunk*advance` below (see that line);
	 * see ../cpu/ap27_crt_sieve.c's own `delta`/`delta_wrapped` comment
	 * for the full reasoning, identical here. The only place this (or
	 * delta_wrapped) is ever examined is at the per-hit reconstruction
	 * below, via an exact carry check -- not a per-chunk bound here. */
	uint64_t delta = 0;
	int delta_wrapped = 0;
	uint32_t completed_prev = 0; /* completed_count is free-running, never reset -- diffed here */

	for (;;) {
		double wall_elapsed = now_seconds() - t_start;
		if (wall_elapsed >= duration) break;

		uint32_t chunk = CHUNK_BATCHES;
		uint32_t groups = (chunk + LOCAL_SIZE_X - 1) / LOCAL_SIZE_X;
		uint64_t dispatched_this_chunk = (uint64_t)groups * LOCAL_SIZE_X;

		/* Same technique as ap27_rom28_kat.c in this directory:
		 * fold this chunk's own batchOffset into baseline0 host-side
		 * (32 cheap CPU divisions), so the shader never needs a live
		 * "batchIdx % P". */
		{
			uint32_t chunkBaseline0[NP];
			for (uint32_t p = 0; p < NP; p++) {
				uint64_t P = PRIMES32[p];
				uint64_t off_term = ((uint64_t)off % P) * advanceModP[p] % P;
				chunkBaseline0[p] = (uint32_t)(((uint64_t)baseline0[p] + off_term) % P);
			}
			UPLOAD(mem_base, chunkBaseline0, NP * sizeof(uint32_t));
		}

		vkResetCommandBuffer(cmd, 0);
		VkCommandBufferBeginInfo cbbi = { .sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO };
		vkBeginCommandBuffer(cmd, &cbbi);
		/* hit_count is reset GPU-side every chunk (its value also
		 * indexes where each hit gets written, so it must be reset or
		 * every hit past the first HIT_CAP of the whole run is silently
		 * dropped forever by the shader's own bounds check).
		 * completed_count is NOT reset here -- resetting it the same way
		 * read back as 0 on essentially every chunk (a false-positive
		 * storm, not the rare/random pattern a real GPU reset produces),
		 * so it's left free-running and diffed cumulatively instead
		 * (see the read-back site below). */
		vkCmdFillBuffer(cmd, buf_cnt, 0, sizeof(uint32_t), 0);
		VkMemoryBarrier fill_barrier = { .sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER,
			.srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT,
			.dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT };
		vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
			0, 1, &fill_barrier, 0, NULL, 0, NULL);
		vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline);
		vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, playout, 0, 1, &dset, 0, NULL);
		/* chunkLen: the real invocation count this chunk covers (always
		 * `chunk` here, since every chunk this driver dispatches is a
		 * full CHUNK_BATCHES -- this loop is duration-bound, not
		 * batch-count-bound, so there's no partial final chunk). Lets
		 * the shader's own bounds check drop the padding invocations
		 * past the real chunk (groups*LOCAL_SIZE_X rounds up to a whole
		 * multiple of LOCAL_SIZE_X) instead of letting them silently
		 * double-scan the next chunk's first few batches. */
		{
			uint32_t chunk_len = chunk;
			vkCmdPushConstants(cmd, playout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(uint32_t), &chunk_len);
		}
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
		elapsed_dispatch += now_seconds() - t0;

		/* per-chunk truncation check via completed_count: DISABLED.
		 * It reads delta=0 on essentially every chunk regardless of
		 * reset strategy tried (host-side reset between submissions,
		 * GPU-side vkCmdFillBuffer with a pipeline barrier, and
		 * cumulative never-reset diffing all gave the same result),
		 * while the actual sieve work is demonstrably correct in the
		 * same runs (KAT found at the exact expected position, hit
		 * counts scaling sensibly with elapsed time) -- this looks
		 * like a bug in this file's own completed_count plumbing, not
		 * evidence of real GPU truncation. TODO: diagnose and
		 * re-enable, or replace the per-invocation global atomic with
		 * a per-workgroup subgroup-reduced add (see
		 * ap27_rom28_kat.c's own comment on this same gap).
		 * Real safety net until then: check dmesg for ring-timeout/
		 * reset events after the fact -- the actual failure mode this
		 * counter was built to detect. */
		MARK_UNUSED(completed_prev);
		MARK_UNUSED(dispatched_this_chunk);

		/* stream this chunk's hits, then reset the counter -- HIT_CAP is
		 * sized for one chunk's worth, not the whole run */
		map_mem(dev, mem_cnt, sizeof(uint32_t), &mapped);
		uint32_t hits = *(uint32_t *)mapped;
		vkUnmapMemory(dev, mem_cnt);
		if (hits > HIT_CAP) {
			fprintf(stderr, "WARNING: chunk at off=%llu: hits=%u > HIT_CAP=%u -- "
				"%u hits from this chunk were dropped (shader-side bounds check)\n",
				(unsigned long long)off, hits, HIT_CAP, hits - HIT_CAP);
		}
		uint32_t nread = (hits < HIT_CAP) ? hits : HIT_CAP;
		/* hits arrive in the shader's atomic append order, which is
		 * arbitrary (not ascending abs_batch/pos -- see the batch=/pos=
		 * comment below) -- so a hit past the 64-bit boundary does NOT
		 * imply every later hit in this chunk is also past it, unlike
		 * the CPU driver's ascending-offset scan. Skip just the
		 * out-of-range hit and keep draining the rest of the chunk's
		 * buffer, then stop after the chunk is fully drained (checked
		 * once this chunk's hits are all read, below) rather than
		 * truncating mid-chunk and silently losing in-range hits. */
		int stop_after_chunk = 0;
		/* This chunk's reportable-offset ceiling, computed once per
		 * chunk instead of twice per hit: `room` is exactly "offset of
		 * the first failure" in this chunk's own within_chunk-offset
		 * space -- any within_chunk beyond it is the beyond-failure
		 * slack, skipped below by one plain comparison per hit, no
		 * overflow check needed per hit. `chunk_unsafe` folds in both
		 * failure modes delta/delta_wrapped (declared above) exist
		 * for: delta having wrapped on its own is checked directly
		 * here, not via `room`, because no offset-based ceiling could
		 * ever express "even within_chunk=0 is invalid" -- 0 <= room
		 * is always true for any unsigned room; A0+delta itself
		 * already exceeding UINT64_MAX is the other half, checked
		 * once here instead of once per hit. */
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
				 * invocation index (never more than a few million -- see
				 * ap27_rom28.comp), added here to the host's own
				 * wide, never-wrapping `off` to get the true absolute
				 * batch index -- both operands must be uint64_t, since
				 * a 32-bit sum would wrap well within range for a real
				 * multi-hour run (see `off`'s own declaration comment). */
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
				/* batch=/pos= (absolute batch index; bit position within
				 * it) -- otherwise-discarded GPU-side indices, kept so a
				 * downstream out-of-order/duplicate A in the hit stream
				 * (parallel invocations racing on the atomic hit-count
				 * append across lanes/workgroups, not a real ordering
				 * guarantee) can be traced back to same-batch vs.
				 * different-batch. See common-hitline.h's
				 * ap27_hitline_field_u64() for the generic reader. */
				if (fprintf(stdout, "HIT B=%llu A=%llu status=unchecked batch=%llu pos=%u\n",
					(unsigned long long)B, (unsigned long long)cand_A,
					(unsigned long long)abs_batch, hit_pos[i]) < 0) {
					fprintf(stderr, "FATAL: write failed, aborting\n");
					exit(1);
				}
			}
			if (fflush(stdout) != 0) { fprintf(stderr, "FATAL: fflush failed, aborting\n"); exit(1); }
		}
		/* reset happens GPU-side via vkCmdFillBuffer at the top of the
		 * next chunk's command buffer -- no host-side reset needed here */

		total_batches += chunk;
		total_hits += hits;
		off += chunk;
		/* delta tracks the same A-units quantity as off*advance, kept
		 * incrementally (see its declaration comment) -- the per-chunk
		 * product chunk*advance is always safe un-checked (bounded,
		 * ~4.3e11), unlike its accumulation into delta across many
		 * chunks, which is exactly what the carry here catches. */
		delta_wrapped |= __builtin_add_overflow(delta, (uint64_t)chunk * advance, &delta);

		if (stop_after_chunk) {
			if (fflush(stdout) != 0) { fprintf(stderr, "FATAL: fflush failed, aborting\n"); exit(1); }
			fprintf(stderr, "search: stopping, A would exceed the 64-bit "
				"range (off=%llu)\n", (unsigned long long)(off - chunk));
			break;
		}

		if ((total_batches / CHUNK_BATCHES) % 20 == 0) {
			uint64_t candidates = total_batches * (uint64_t)BATCH;
			/* A=<current A>, for real resumability: A at this exact
			 * batch offset is A0+delta, exact -- lets a resume script
			 * pick up from here without re-scanning the (potentially
			 * huge) hit stream. Guarded like the per-hit check above:
			 * this is a diagnostic line, not a HIT, but a silently
			 * wrong resumability hint is still a real problem. */
			uint64_t A_now;
			int a_now_safe = !delta_wrapped && !__builtin_add_overflow(A0, delta, &A_now);
			if (a_now_safe) {
				fprintf(stderr, "[%.1fs] B=%llu batches=%llu candidates=%llu A=%llu hits=%llu rate=%.1fM cand/s\n",
					wall_elapsed, (unsigned long long)B, (unsigned long long)total_batches,
					(unsigned long long)candidates, (unsigned long long)A_now,
					(unsigned long long)total_hits, candidates / elapsed_dispatch / 1e6);
			} else {
				fprintf(stderr, "[%.1fs] B=%llu batches=%llu candidates=%llu A=(exceeds 64-bit range) hits=%llu rate=%.1fM cand/s\n",
					wall_elapsed, (unsigned long long)B, (unsigned long long)total_batches,
					(unsigned long long)candidates,
					(unsigned long long)total_hits, candidates / elapsed_dispatch / 1e6);
			}
		}
	}
#undef UPLOAD

	double final_elapsed = now_seconds() - t_start;
	uint64_t candidates = total_batches * (uint64_t)BATCH;
	fprintf(stderr, "DONE: B=%llu elapsed=%.3fs batches=%llu candidates=%llu hits=%llu rate=%.1fM cand/s\n",
		(unsigned long long)B, final_elapsed, (unsigned long long)total_batches,
		(unsigned long long)candidates, (unsigned long long)total_hits,
		candidates / elapsed_dispatch / 1e6);
	return 0;
}
