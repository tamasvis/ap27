/* ap27_crt_sieve.c -- CRT/ROM-batch sieve for AP27 candidates, CPU reference.
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * This is the minimal, standalone implementation of the design: no GPU,
 * no external dependencies beyond libc. If you want to understand the
 * mechanism before looking at the GPU port (../gpu/), start here. See
 * the top-level README for the full design rationale (why residue-
 * indexed bitmasks beat recomputed ones, why 32 primes, the small-prime
 * CRT-merge optimization used below) -- this comment only summarizes
 * enough to read the code.
 *
 * A candidate's hot-path representation is its residue tuple mod each
 * of a fixed set of primes (a CRT representation, not an integer),
 * reconstructed to a real A only once it survives every prime check.
 * Candidate positions follow a fixed period of 5394 (= 6*29*31) per B,
 * so the next 512 candidate offsets from any starting point are known
 * in advance -- meaning "is candidate s (of the next 512) divisible by
 * P" is a fixed 512-bit answer that depends ONLY on the current residue
 * mod P, never on absolute position. Precompute that mask for every
 * residue (one table of size P per prime); the per-batch sieve step is
 * then just one lookup per prime, AND them all together -- no division
 * or modulo anywhere in the hot loop. Advancing to the next batch is a
 * single modular add per prime, not a recomputation of anything.
 *
 * Below: the ROM tables for the 32-prime design (7 primes dividing C,
 * CRT-merged into 3 combined-modulus groups -- 28 lookups, not 32 --
 * plus the 25 table primes 37..151), built once per B, then a
 * sieve loop that scans forward from a known AP27 (a KAT, known-
 * answer-test) and confirms it's found at the exact expected position
 * -- the standard correctness check this whole project uses throughout.
 * See KAT_A[]/KAT_B[] below for the real, publicly verifiable examples
 * used.
 *
 * Build:   make (from the repo root) -- see the Makefile for the exact
 *          flags, rather than duplicating them here where they could
 *          drift out of sync.
 * Run:     ./ap27_crt_sieve [kat_index 1..N] [window]
 *          ./ap27_crt_sieve --timed B=<b> A=<a_lo> SECONDS=<duration>
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include "../common-modinv.h"
#include "../common-cliargs.h"

#define BATCH   512u
#define WORDS   (BATCH / 32u)   /* 16 uint32_t words per 512-bit mask */
#define NB      25u             /* the 25 primes 37..151 that don't divide C */
#define NGROUP  3u              /* CRT-merged small-prime groups */
#define GROUP_MAX 3u            /* largest group has 3 primes */

/* The 25 primes 37..151 that don't divide C -- each still needs its own
 * real per-position lookup. (A candidate's residue mod P plus the fixed
 * step BC mod P does determine all 27 terms mod P; the real distinction
 * from 5..23 is that for P not dividing C, which of the 27 terms hits P
 * actually varies by position, so there's something to look up at all --
 * for P | C every term is always == A, collapsing to one check. Merging
 * these 25 the same way 5..23 are merged was considered and rejected:
 * table growth scales with the *product* of grouped primes, too costly
 * once the primes themselves aren't tiny.)
 *
 * TABLE_PRIME_ADVANCE[] below steps each one's residue index forward
 * by one batch using a fixed integer increment plus a conditional
 * wraparound (the same idea as Bresenham's line algorithm) -- never a
 * real division.
 */
static const uint16_t TABLE_PRIMES[NB] = {
	37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
	101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151
};
/* advance % TABLE_PRIMES[i]; advance = 32 * 6 * 29 * 31
 * each bitmask advances 32*16 mod-6*29*31 full cycles
 */
static const uint16_t TABLE_PRIME_ADVANCE[NB] = {
        3, 39, 6, 24, 40, 33, 39, 16, 7, 36, 72, 51, 37, 45,
        100, 83, 17, 61, 57, 15, 81, 125, 109, 66, 15
};

/* "Tough luck" B rejection: B sharing a factor with a small prime this
 * design relies on breaks table construction, via two genuinely
 * DIFFERENT mechanisms:
 *
 * (1) B % <table prime, 37..151> == 0: the closed-form
 * D = (P - (B%P*C%P)%P) % P used to build that prime's table is IDENTICALLY
 * ZERO -- same family as the well-known "P|C" degeneracy, just triggered by
 * B instead of by P. With D=0, "good = ((raw*Dinv)%P) > 26" would be
 * "0 > 26" for every position -- that one prime's table permanently
 * "always bad", ANDed into every batch, silently rejecting 100% of
 * candidates for the rest of the run, with no throughput anomaly to
 * notice it by. Caught here, before it can happen: a known-bad B is
 * rejected with a clear message below, rather than reaching
 * modinv_u64(0,P) at all. (If it ever did, that call's own
 * post-condition check would abort too -- see common-modinv.h -- not
 * silently produce this D=0 case; this guard exists for a clearer,
 * earlier error, not because modinv_u64 itself would get it wrong.)
 *
 * (2) B % 29 == 0 or B % 31 == 0: a DIFFERENT, WORSE failure mode.
 * build_good5394()'s g29/g31 arrays encode "the legal starting
 * residues mod 29/31" as g29={step29, 2*step29 mod 29} where
 * step29=(B*C) mod 29 -- valid only when step29!=0 (B%29!=0), since then
 * the 27 terms occupy 27 DISTINCT residues mod 29, leaving exactly 2 safe
 * ones. When B%29==0, step29=0, so ALL 27 terms collapse to the SAME
 * residue (A mod 29) -- meaning the real safe set is 28 of 29 residues
 * (everything except A%29==0), the total opposite of what g29={0,0} (the
 * formula's actual output when step29=0) tells the gate to search for.
 * The gate then ACTIVELY GENERATES ONLY A%29==0 positions -- the one
 * genuinely bad residue -- as if they were legal. Since 29 (and 31) are
 * NOT among the 32 primes this design's own sieve checks (both are
 * handled ENTIRELY by this gate, by design), every resulting "survivor"
 * sails through the 32-prime sieve undetected -- not reduced coverage,
 * not zero hits, but a 100% silent false-positive rate: every single
 * reported hit is trivially divisible by 29 or 31, while throughput and
 * hit counts both look completely normal. Arguably more dangerous than
 * case (1) precisely because it does NOT look anomalous.
 *
 * ~32.8% of uniformly-chosen B's fail this check overall -- refusing
 * loudly at init, rather than building real handling into either closed
 * form, is the deliberate choice here. */
static int ap27_b_is_tough_luck(uint64_t B, const char **why_prime, const char **why_mechanism)
{
	if (B % 29 == 0) {
		if (why_prime) *why_prime = "29";
		if (why_mechanism) *why_mechanism = "good5394 gate";
		return 1;
	}
	if (B % 31 == 0) {
		if (why_prime) *why_prime = "31";
		if (why_mechanism) *why_mechanism = "good5394 gate";
		return 1;
	}
	for (uint32_t i = 0; i < NB; i++) {
		if (B % TABLE_PRIMES[i] == 0) {
			static char buf[16];
			snprintf(buf, sizeof(buf), "%u", TABLE_PRIMES[i]);
			if (why_prime) *why_prime = buf;
			if (why_mechanism) *why_mechanism = "table-prime Dinv";
			return 1;
		}
	}
	return 0;
}

static void ap27_reject_tough_luck_b(uint64_t B)
{
	const char *why_prime = NULL, *why_mechanism = NULL;
	if (ap27_b_is_tough_luck(B, &why_prime, &why_mechanism)) {
		if (strcmp(why_mechanism, "good5394 gate") == 0) {
			fprintf(stderr,
				"FAIL: B=%llu is divisible by %s -- this gate's own closed form "
				"(build_good5394) degenerates when B shares this factor: it would "
				"silently generate ONLY A%%%s==0 starting positions (the one "
				"genuinely composite residue, not a legal one), and since %s isn't "
				"checked by the 32-prime sieve itself, every resulting \"survivor\" "
				"would be a silent, 100%% false positive with completely normal-"
				"looking throughput. Not implemented as a special case "
				"(\"tough luck\") -- choose a B not divisible by 29, 31, or "
				"37..151.\n",
				(unsigned long long)B, why_prime, why_prime, why_prime);
		} else {
			fprintf(stderr,
				"FAIL: B=%llu is divisible by the table prime %s -- D=0 in "
				"that prime's table-build closed form, which would make that "
				"prime's table permanently \"always bad\" -- zero real hits "
				"for the whole run, despite completely normal-looking "
				"throughput. Not implemented as a special case "
				"(\"tough luck\") -- choose a B not divisible by 29, 31, "
				"or 37..151.\n",
				(unsigned long long)B, why_prime);
		}
		exit(1);
	}
}

/* The 7 primes dividing C, CRT-merged into 3 groups (one triple, two
 * pairs): {5,7,11}, {13,23}, {17,19}. GROUP_LEN[k] entries of GROUP[k][]
 * are used; the rest are padding (value 1, a no-op factor) so every row
 * has the same width. */
static const uint32_t GROUP[NGROUP][GROUP_MAX] = { {5,7,11}, {13,23,1}, {17,19,1} };
static const uint32_t GROUP_PRODUCTS[NGROUP] = { 385, 299, 323 };
static const uint32_t GROUP_LEN[NGROUP] = { 3, 2, 2 };
/**/
#define MAX_GROUP_ELEMS        ((unsigned int) 5 * 7 * 11)
/* nr. elems with padding forcing to multiple of 64 */
#define MAX_GROUP_ELEMS_PADDED ((unsigned int) 448)


static const uint64_t C_VAL = 223092870ull;

/* Real, publicly verifiable known AP27 examples. Each is a genuine
 * 27-term arithmetic progression of primes; the sieve must find each one
 * exactly, at the exact predicted position, or something is wrong.
 */
static const uint64_t KAT_A[] = {
	605185576317848261ull,
	277699295941594831ull,
	224584605939537911ull,
	973044677408627339ull,
};
static const uint64_t KAT_B[] = {
	155368778ull,
	170826477ull,
	81292139ull,
	174505864ull,
};
#if !defined(ARRAY_ELEMS)
#define  ARRAY_ELEMS(arr)  (sizeof(arr) / sizeof((arr)[0]))
#endif
#define NKATS ((int)ARRAY_ELEMS(KAT_A))

/* The mod-5394 (= 6*29*31) legal-position cycle for a given B: which of
 * the 5394 residue classes mod 5394 are even candidates for a legal
 * (6k+1)/(6k+5) start whose B-mod-29/31 class doesn't get excluded
 * outright. This cycle is fixed for the life of a search on this B (16
 * legal positions per 5394-wide block, summing to exactly 5394), which is
 * what makes the whole batch-of-512 affine offset trick below possible.
 * uint32_t rather than the uint16_t these values (all < 5394) would fit
 * in: this runs during initialization only, so the few bytes saved
 * aren't worth a narrower type that C's integer promotion rules would
 * widen back to int in every expression anyway. */
static uint32_t good5394[16];
static int gap5394[16];

static void build_good5394(uint64_t b)
{
	uint64_t bc29 = (b % 29) * C_VAL % 29;
	uint64_t bc31 = (b % 31) * C_VAL % 31;
	uint8_t g29[2] = { (uint8_t)bc29, (uint8_t)((2 * bc29) % 29) };
	uint8_t g31[4] = { (uint8_t)bc31, (uint8_t)((2 * bc31) % 31),
	                    (uint8_t)((3 * bc31) % 31), (uint8_t)((4 * bc31) % 31) };
	uint32_t combined[16];
	int n = 0;
	int mod6[2] = { 1, 5 };
	for (int mi = 0; mi < 2; mi++) {
		for (int gi = 0; gi < 2; gi++) {
			for (int hi = 0; hi < 4; hi++) {
				uint32_t r;
				for (r = 0; r < 5394; r++)
					if (r % 6 == (unsigned)mod6[mi] && r % 29 == g29[gi] && r % 31 == g31[hi])
						break;
				combined[n++] = r;
			}
		}
	}
	/* insertion sort into ascending order */
	for (int i = 1; i < 16; i++) {
		uint32_t key = combined[i];
		int j = i - 1;
		while (j >= 0 && combined[j] > key) { combined[j + 1] = combined[j]; j--; }
		combined[j + 1] = key;
	}
	for (int i = 0; i < 16; i++) good5394[i] = combined[i];
	for (int i = 0; i < 16; i++) gap5394[i] = (int)((combined[(i + 1) % 16] - combined[i] + 5394) % 5394);
}

/* rel_offset[s] = the offset from the batch's own starting A to the s-th
 * of the next 512 legal candidate positions, for s in [0, 512). Built
 * once per B from the 16-entry gap cycle above (32 repeats of the cycle
 * = 512 positions), independent of where in A-space the batch actually
 * starts -- this is the affine relationship that makes the whole ROM
 * table position-independent. */
static uint64_t rel_offset[BATCH];

static void build_rel_offsets(void)
{
	uint64_t cum = 0;
	int idx = 0;
	for (int block = 0; block < 32; block++)
		for (int g = 0; g < 16; g++) { rel_offset[idx++] = cum; cum += (uint64_t)gap5394[g]; }
}

/* Two ROM table groups, one per mechanism:
 *  - rom_group[k]: the CRT-merged small-prime groups, indexed mod the
 *    group's product, a good bit means none of the group's primes divide
 *    this position.
 *  - rom_b[p]: the 25 real per-position tables, one per table prime
 * Both are pure functions of residue, built once per B, never touched
 * again. */
static uint32_t rom_group[NGROUP][MAX_GROUP_ELEMS_PADDED*WORDS];
static uint32_t *rom_b[NB];
static uint64_t Dinv_b[NB];

static void build_rom(uint64_t B)
{
	for (uint32_t k = 0; k < NGROUP; k++) {
		uint32_t prod = GROUP_PRODUCTS[k];
		for (uint32_t base = 0; base < prod; base++) {
			uint32_t *mask = &rom_group[k][base * WORDS];
			for (uint32_t s = 0; s < BATCH; s++) {
				uint64_t raw = (base + rel_offset[s]) % prod;
				/* good iff A is divisible by NONE of this group's primes --
				 * the OR of the group's individual "bad" checks, baked
				 * into one table. */
				int good = 1;
				for (uint32_t i = 0; i < GROUP_LEN[k]; i++)
					if (raw % GROUP[k][i] == 0) good = 0;
				if (good) mask[s / 32] |= (1u << (s % 32));
			}
		}
	}

	for (uint32_t p = 0; p < NB; p++) {
		uint16_t P = TABLE_PRIMES[p];
		uint16_t c = C_VAL % P;
		uint64_t r1 = B % P;
		uint64_t D = (P - (r1 * c) % P) % P;
		uint64_t D_inv = modinv_u64(D, P);
		Dinv_b[p] = D_inv;
		rom_b[p] = calloc((size_t)P * WORDS, sizeof(uint32_t));
		for (uint64_t base = 0; base < P; base++) {
			uint32_t *mask = &rom_b[p][base * WORDS];
			for (uint32_t s = 0; s < BATCH; s++) {
				uint64_t raw = (base + rel_offset[s]) % P;
				/* "term index i satisfies P | term" reduces to a single
				 * closed-form comparison via the precomputed inverse --
				 * the key property used here is just that it's a pure
				 * function of raw. */
				int good = ((raw * D_inv) % P) > 26;
				if (good) mask[s / 32] |= (1u << (s % 32));
			}
		}
	}
}

// allow SIMDification without considerations for loop unrolling etc.
static inline void and_words_inplace(uint32_t dst[restrict static WORDS],
                               const uint32_t src[restrict static WORDS])
{
	_Static_assert(WORDS == 16, "and_words_inplace is manually unrolled for exactly 16 words");
	dst[  0 ] &= src[  0 ];
	dst[  1 ] &= src[  1 ];
	dst[  2 ] &= src[  2 ];
	dst[  3 ] &= src[  3 ];
	dst[  4 ] &= src[  4 ];
	dst[  5 ] &= src[  5 ];
	dst[  6 ] &= src[  6 ];
	dst[  7 ] &= src[  7 ];
	dst[  8 ] &= src[  8 ];
	dst[  9 ] &= src[  9 ];
	dst[ 10 ] &= src[ 10 ];
	dst[ 11 ] &= src[ 11 ];
	dst[ 12 ] &= src[ 12 ];
	dst[ 13 ] &= src[ 13 ];
	dst[ 14 ] &= src[ 14 ];
	dst[ 15 ] &= src[ 15 ];
}

static inline uint32_t u32_ge_or0(uint32_t val, uint32_t threshold)
{
        return ((val >= threshold) ? threshold : 0);
}

/*
 * run_timed_search(duration_s, B, A_start) -- generic, open-ended
 * search: blindly walks forward from A_start along fixed B for
 * duration_s seconds, printing every sieve-surviving candidate.
 *
 * Compact line shape: "HIT B=..<5-digit suffix of B> A=x<hex, no
 * leading zero> status=unchecked batch=... pos=... id=<5-digit stream
 * id>". Full B is logged once already, in the "ap27_crt_sieve --timed:
 * B=..." header line printed at the top of this function -- every
 * subsequent HIT line then only needs a short suffix for eyeball-
 * disambiguation of interleaved concurrent-stream logs (see
 * common-hitline.h's ap27_print_b_suffix/ap27_parse_u128_auto/
 * ap27_hitline_parse_compact for the shared grammar/parsing side of
 * this). The GPU driver (../gpu/ap27_rom28_search.c) still emits
 * full decimal "B=... A=..." per line, so the two genuinely diverge in
 * format. Surviving this 32-prime sieve is not confirmation, just a
 * candidate for a downstream second-stage filter (presieve/) and,
 * eventually, a real Fermat/Miller-Rabin test (not included here).
 * batch=/pos= are this design's own indices (batch = 512-candidate
 * ROM lookup round, pos = bit position within it), same diagnostic
 * role as the GPU driver's fields, though a single-threaded sequential
 * CPU scan has no reason to ever produce them out of order the way
 * parallel GPU lanes can. id= is this process's own PID mod 1e5 --
 * distinct concurrent workers (e.g. several instances of this binary
 * launched against different A0's for the same B) can share one B
 * while scanning entirely different A ranges, so downstream per-stream
 * gap-tracking needs a key that doesn't assume B alone identifies one
 * worker's stream.
 *
 * CLI contract: takes B/A0/duration as dd-style "B=.. A=.. SECONDS=.."
 * tokens, in any order (see ../common-cliargs.h), same keys as the GPU
 * search driver. A caller can launch several instances (one process
 * per B/A0 range) with no wrapper logic needed on either side. This
 * driver only ever walks forward (a
 * mirrored backward walk would need subtracting 'advance' mod each
 * modulus instead of adding, not implemented since nothing currently
 * needs it).
 *
 * A_start need not be a legal position: snapped up to the next A with
 * A mod 5394 == good5394[0] specifically (not the nearest of the 16
 * legal residues per period -- this can skip up to 15 closer legal
 * positions within the same period, harmless since the scan is
 * open-ended anyway), same requirement as main()'s own A0 construction
 * below.
 */
static void run_timed_search(uint64_t duration_s, uint64_t B, uint64_t A_start)
{
	build_good5394(B);
	build_rel_offsets();
	build_rom(B);

	uint64_t anchor = (uint64_t)good5394[0];
	uint64_t rem = A_start % 5394ull;
	uint64_t A0 = (rem <= anchor) ? (A_start - rem + anchor)
	                              : (A_start - rem + 5394ull + anchor);

	uint32_t base_g[NGROUP], base_b[NB];
	for (uint32_t k = 0; k < NGROUP; k++) base_g[k] = (uint32_t)(A0 % GROUP_PRODUCTS[k]);
	for (uint32_t p = 0; p < NB; p++) base_b[p] = (uint32_t)(A0 % TABLE_PRIMES[p]);

	/* Tracks the batch offset from A0, not the absolute A -- deliberately
	 * never range-checked on every `+= advance` below; the only place
	 * either of these is ever examined is right where a candidate is
	 * about to be reported, via an exact carry check (not a bound here).
	 * delta_wrapped is a sticky flag (via `|=`, every batch, from that
	 * add's own native carry bit -- not a separate comparison): once
	 * set by any wrap of delta itself, however many, it stays set, so
	 * the two checks together are exact, not probabilistic. */
	uint64_t delta = 0;
	int delta_wrapped = 0;
	uint64_t batches = 0, hits = 0;
	const uint64_t advance = 32ull * 5394ull;

	/* Per-invocation stream id, for downstream gap-tracking to key on:
	 * distinct concurrent workers can share the same B while scanning
	 * entirely different A ranges (e.g. several instances launched
	 * against different A0's for one B), so B alone isn't enough to
	 * tell two streams apart. getpid() is a cheap, already-unique-
	 * among-concurrent-processes choice -- no RNG/collision handling
	 * needed, unlike a hand-rolled random token. The FULL pid, not a
	 * truncated form: a 5-digit-style modulus (100,000) would let two
	 * real, concurrently-live PIDs that happen to be 100,000 apart
	 * collide -- easily possible, since pid_max is often 4,194,304 on
	 * Linux -- exactly defeating the point of using getpid() at all.
	 * The id= field parses generically by key (see
	 * ap27_hitline_field_u64()), so there is no fixed-width wire
	 * format to preserve here. */
	unsigned stream_id = (unsigned)getpid();

	time_t t0 = time(NULL), tend = t0 + (time_t)duration_s;
	printf("ap27_crt_sieve --timed: B=%llu A0=%llu duration=%llus id=%u\n",
		(unsigned long long)B, (unsigned long long)A0, (unsigned long long)duration_s, stream_id);
	fflush(stdout);

	for (;;) {
		uint32_t combined[WORDS];

		memcpy(combined, &(rom_group[0][base_g[0] * WORDS]), sizeof(combined));
		and_words_inplace(combined, &(rom_group[1][base_g[1] * WORDS]));
		and_words_inplace(combined, &(rom_group[2][base_g[2] * WORDS]));
		for (uint32_t p = 0; p < NB; p++)
			and_words_inplace(combined, &(rom_b[p][base_b[p] * WORDS]));

		for (int w = 0; w < (int)WORDS; w++) {
			uint32_t bits = combined[w];
			while (bits != 0) {
				int b = __builtin_ctz(bits);
				int s = w * 32 + b;
				/* The only overflow check in this whole function, run
				 * only for an actual candidate -- not every batch. Two
				 * independent layers, together exact (not
				 * probabilistic): `delta_wrapped` is a sticky
				 * carry flag from `delta`'s own `+= advance` below (see
				 * that line -- branch-free, an OR of the add's native
				 * carry bit, so it costs nothing extra there), catching
				 * the case where delta itself has wrapped past 2^64 at
				 * some earlier batch with no hit in between to report
				 * it; the add_overflow chain catches the much more
				 * reachable case where A0+delta already exceeds
				 * UINT64_MAX while delta has not yet wrapped on its own
				 * (always possible first, since A0 can be up to
				 * UINT64_MAX/2). Either one stops the run here, cleanly,
				 * instead of printing a wrapped A. */
				uint64_t cand_A, t;
				if (delta_wrapped ||
					__builtin_add_overflow(delta, rel_offset[s], &t) ||
					__builtin_add_overflow(A0, t, &cand_A)) {
					fflush(stdout);
					fprintf(stderr, "ap27_crt_sieve --timed: stopping, A would "
						"exceed the 64-bit range at batch=%llu\n",
						(unsigned long long)batches);
					goto done;
				}
				/* format fully into a stack buffer, THEN one fwrite +
				 * fflush -- N instances of this binary can run
				 * concurrently, all sharing the SAME downstream pipe.
				 * Without flushing every line, libc's
				 * default full-buffering-on-a-pipe batches many lines
				 * into one write(), and once that exceeds PIPE_BUF
				 * (4096B on Linux) the kernel no longer guarantees
				 * atomicity -- another process's write can interleave
				 * mid-flush, tearing a line in two (real, observed
				 * corruption, not just reordering: "status=unc221
				 * status=unchecked"). A latent race, visible only once
				 * per-worker throughput is high enough to fill a
				 * multi-KB buffer between flushes -- which it routinely
				 * is here. Building the line in memory first (snprintf)
				 * then writing it whole in one fwrite+fflush keeps every
				 * actual write() a single, small (~70-90 byte),
				 * PIPE_BUF-safe atomic unit, regardless of thread count
				 * or throughput. A torn line's B/A would still be
				 * recoverable by re-scanning that exact range, but
				 * writing it correctly up front avoids needing to. */
				{
					char linebuf[128];
					int n = snprintf(linebuf, sizeof(linebuf),
						"HIT B=..%05llu A=x%llx status=unchecked batch=%llu pos=%d id=%u\n",
						(unsigned long long)(B % 100000ull), (unsigned long long)cand_A,
						(unsigned long long)batches, s, stream_id);
					fwrite(linebuf, 1, (size_t)n, stdout);
					fflush(stdout);
				}
				hits++;
				bits &= bits - 1;
			}
		}

		/* literal constants, not GROUP_PRODUCTS[k] -- required for the
		 * compiler to see a compile-time-constant modulus and apply
		 * Granlund-Montgomery strength reduction (an array read
		 * defeats it, or at least isn't guaranteed not to); hardwired
		 * the same way main()'s own specialized loop already is
		 * below, regenerate both if the grouping ever changes. */
		base_g[0] += advance; base_g[0] %= 385;   /* 5 * 7 * 11 */
		base_g[1] += advance; base_g[1] %= 299;   /* 13 * 23 */
		base_g[2] += advance; base_g[2] %= 323;   /* 17 * 19 */
		for (uint32_t p = 0; p < NB; p++) {
			base_b[p] += TABLE_PRIME_ADVANCE[p];
			base_b[p] -= u32_ge_or0(base_b[p], TABLE_PRIMES[p]);
		}
		/* Unconditional, unchecked, every batch -- the entire point of
		 * deferring the real check to the hit site above. The OR into
		 * delta_wrapped is the add's own native carry bit, not a
		 * separate comparison -- no branch added to this path. */
		delta_wrapped |= __builtin_add_overflow(delta, advance, &delta);
		batches++;

		if ((batches & 0xFFFFFull) == 0) {
			time_t now = time(NULL);
			double elapsed = difftime(now, t0);
			double rate = (elapsed > 0) ? (double)(batches * BATCH) / elapsed : 0.0;
			printf("PROGRESS t=%llds batches=%llu candidates=%llu hits=%llu rate=%.0f/s\n",
				(long long)elapsed, (unsigned long long)batches,
				(unsigned long long)(batches * BATCH),
				(unsigned long long)hits, rate);
			fflush(stdout);
			if (now >= tend) break;
		}
	}

done:
	printf("TIMED_DONE: B=%llu batches=%llu candidates=%llu hits=%llu\n",
		(unsigned long long)B, (unsigned long long)batches,
		(unsigned long long)(batches * BATCH), (unsigned long long)hits);
	fflush(stdout);

	for (uint32_t p = 0; p < NB; p++)
		free(rom_b[p]);
}

int main(int argc, char **argv)
{
	if (argc > 1 && strcmp(argv[1], "--report-primes") == 0) {
		printf("32\n");
		return 0;
	}
	if (argc > 1 && strcmp(argv[1], "--timed") == 0) {
		const char *usage =
			"usage: ap27_crt_sieve --timed B=<b> A=<a0> SECONDS=<secs> (any order)\n";
		uint64_t B = 0, A_start = 0, duration_s = 0;
		ap27_kv_spec specs[] = {
			/* A's ceiling (UINT64_MAX/2) is not a representability limit
			 * (see run_timed_search's own comment on the 64-bit stop) --
			 * it's generously far below 2^64 so a bad/adversarial A_lo
			 * fails here, immediately and by name, rather than running
			 * first and hitting that stop microseconds later. */
			{ .key = "B",       .type = AP27_KV_U64, .u64_out = &B,         .min = 1, .max = UINT64_MAX,   .required = 1 },
			{ .key = "A",       .type = AP27_KV_U64, .u64_out = &A_start,   .min = 0, .max = UINT64_MAX / 2, .required = 1 },
			{ .key = "SECONDS", .type = AP27_KV_U64, .u64_out = &duration_s, .min = 0, .max = 100000000ull, .required = 1 },
		};
		if (!ap27_parse_kv_args(argc, argv, 2, specs, 3, usage))
			return 1;
		ap27_echo_kv_args("ap27_crt_sieve --timed", specs, 3);
		ap27_reject_tough_luck_b(B);
		run_timed_search(duration_s, B, A_start);
		return 0;
	}

	/* dd-style "KEY=value" args here too, same reasoning as --timed
	 * above, via the shared parser in ../common-cliargs.h. */
	const char *katusage = "usage: ap27_crt_sieve [KAT=<1..4>] [WINDOW=<n>] (any order)\n";
	uint64_t kat_idx64 = 1, window = 3000000000000ull;
	ap27_kv_spec katspecs[] = {
		{ .key = "KAT",    .type = AP27_KV_U64, .u64_out = &kat_idx64, .min = 1, .max = (uint64_t)NKATS },
		{ .key = "WINDOW", .type = AP27_KV_U64, .u64_out = &window,    .min = 1, .max = 10000000000000000000ull },
	};
	if (!ap27_parse_kv_args(argc, argv, 1, katspecs, 2, katusage))
		return 1;
	ap27_echo_kv_args("ap27_crt_sieve", katspecs, 2);
	int kat_idx = (int)kat_idx64;
	uint64_t B = KAT_B[kat_idx - 1];
	uint64_t KAT_TARGET = KAT_A[kat_idx - 1];

	ap27_reject_tough_luck_b(B);
	build_good5394(B);
	build_rel_offsets();
	build_rom(B);

	size_t rom_bytes = 0;
	for (uint32_t k = 0; k < NGROUP; k++) rom_bytes += (size_t)GROUP_PRODUCTS[k] * WORDS * sizeof(uint32_t);
	for (uint32_t p = 0; p < NB; p++) rom_bytes += (size_t)TABLE_PRIMES[p] * WORDS * sizeof(uint32_t);
	printf("ap27_crt_sieve: KAT #%d, B=%llu, target A=%llu\n",
		kat_idx, (unsigned long long)B, (unsigned long long)KAT_TARGET);
	printf("ROM table: %u merged groups + %u table primes = %u lookups/batch, "
		"%zu bytes (%.1f KB)\n", NGROUP, NB, NGROUP + NB, rom_bytes, rom_bytes / 1024.0);

	/* WINDOW has no meaningful static ceiling of its own -- it has to
	 * stay well below whichever KAT's target this run picked, which
	 * isn't known until kat_idx is. Checked here, by name, rather than
	 * letting `KAT_TARGET - window` underflow into a huge k_ below and
	 * surface only as a confusing "KAT NOT FOUND" much later. */
	if (window + (uint64_t)good5394[0] >= KAT_TARGET) {
		fprintf(stderr, "FAIL: WINDOW=%llu is too large for KAT #%d's target "
			"A=%llu (must stay well below it)\n",
			(unsigned long long)window, kat_idx, (unsigned long long)KAT_TARGET);
		return 1;
	}

	/* Construct A0 directly as a legal position -- good5394[0] plus
	 * enough whole 5394-cycles to land just below the target -- rather
	 * than starting from an arbitrary integer and correcting afterward.
	 * A0 mod 5394 == good5394[0] is required for rel_offset[] to be
	 * valid; this is what guarantees it by construction. */
	uint64_t k_ = (KAT_TARGET - window - (uint64_t)good5394[0]) / 5394ull;
	uint64_t A0 = (uint64_t)good5394[0] + k_ * 5394ull;

	uint32_t base_g[NGROUP], base_b[NB];
	for (uint32_t k = 0; k < NGROUP; k++) base_g[k] = (uint32_t)(A0 % GROUP_PRODUCTS[k]);
	for (uint32_t p = 0; p < NB; p++) base_b[p] = (uint32_t)(A0 % TABLE_PRIMES[p]);
		/* note: divisions outside loop are fine */

	uint64_t A = A0;
	uint64_t batches = 0;
	int found = 0;
	const uint64_t advance = 32ull * 5394ull; /* one 512-candidate batch's worth of A */

	while (A < KAT_TARGET + 200000000ull) {
		uint32_t combined[WORDS];

		/* Sequential AND chain, not a hand-rolled binary-tree
		 * reduction: the compiler already autovectorizes this into
		 * wide SIMD ANDs on every target checked (x86-64 AVX-512,
		 * ARM64 NEON -- see the README's SIMD section), so
		 * restructuring by hand would need to beat that with a real
		 * measurement first. */
			// AND(all-1, first bitmask) == first bitmask
		memcpy(combined, &(rom_group[0][base_g[0] * WORDS]),
		       sizeof(combined));

		and_words_inplace(combined, &(rom_group[1][base_g[1] * WORDS]));
		and_words_inplace(combined, &(rom_group[2][base_g[2] * WORDS]));

		for (uint32_t p = 0; p < NB; p++) {
			const uint32_t *m = &rom_b[p][base_b[p] * WORDS];
			and_words_inplace(combined, m);
		}

		for (int w = 0; w < (int)WORDS; w++) {
			uint32_t bits = combined[w];
			while (bits != 0) {
				int b = __builtin_ctz(bits);
				int s = w * 32 + b;
				uint64_t cand_A = A + rel_offset[s];
				if (cand_A == KAT_TARGET) {
					printf("FOUND: A=%llu (batch=%llu, pos=%d) -- matches KAT #%d exactly\n",
						(unsigned long long)cand_A, (unsigned long long)batches, s, kat_idx);
					found = 1;
				}
				bits &= bits - 1;
			}
		}

                        /* base_g[i] += advance   mod GROUP_PRODUCTS[i] */
                        /* directly inlined constants allow */
                        /* Granlund-Montgomery specialization */
                base_g[0] += advance;
                base_g[0] %= 385;   /* 5 * 7 * 11 */
                base_g[1] += advance;
                base_g[1] %= 299;   /* 13 * 23 */
                base_g[2] += advance;
                base_g[2] %= 323;   /* 17 * 19 */

                                /* modular add; no need for real remainder */
                for (uint32_t p = 0; p < NB; p++) {
                        base_b[p] += TABLE_PRIME_ADVANCE[p];
                        base_b[p] -= u32_ge_or0(base_b[p], TABLE_PRIMES[p]);
                }
                A += advance;
                batches++;
        }

	for (uint32_t p = 0; p < NB; p++)
		free(rom_b[p]);

	printf("scanned %llu batches (%llu candidates), KAT %s\n",
		(unsigned long long)batches, (unsigned long long)(batches * BATCH),
		found ? "FOUND" : "NOT FOUND -- something is wrong");

	return !found;
}
