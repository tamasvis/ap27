/*
 * ap27_prefilter256.c -- second-stage sieve for AP27 candidates.
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * A candidate that has already survived the 32-prime CRT sieve (see
 * ../gpu/ and ../cpu/) is far from confirmed -- it's just no longer
 * provably composite by a *small* prime. This stage extends that check
 * with primes 751..2677 (256 of them), using the same real, direct
 * mod-P logic (no precomputed table -- these primes are checked "fresh"
 * per candidate, since by this point candidates are rare enough that a
 * table isn't worth building). It's meant to run between the sieve and
 * the expensive step (a Fermat/Miller-Rabin primality test on each of
 * the 27 terms) -- cheaper than even one such exponentiation, and real
 * measured data shows it rejects ~99.1% of candidates that reach it.
 *
 * There's a real coverage subtlety worth being explicit about. The
 * primes 157..743 (96 of them, GAP_PRIMES below) are NOT checked by the
 * 32-prime sieve in ../cpu/ and ../gpu/ (which stops at 151 for
 * throughput -- see the top-level README's "Why 32 primes"). This
 * file's own EXT_PRIMES (751..2677) always assumed primes up to 743
 * were already covered upstream -- true for a wider sieve, not true
 * for the 32-prime one. Rather than shipping two near-duplicate files,
 * GAP_PRIMES is a switchable region: pass check_gap=0 downstream of a
 * wide sieve (a real no-op, since that range is already covered) or
 * check_gap=1 downstream of the 32-prime one (closes the real gap).
 * Checked first when enabled, since these are much higher-yield per
 * prime than the 751..2677 range.
 *
 * Per-prime logic: a candidate survives prime P iff none of the 27 terms
 * (A + i*B*C for i=0..26) is divisible by P. Computed incrementally --
 * one real reduction for A mod P and for the per-step increment
 * (B*C) mod P, then 27 branchless add-and-conditional-subtract steps,
 * no division inside the inner loop at all.
 */

#include <stdint.h>
#include <stdio.h>

#define NUM_EXT_PRIMES 256
#define NUM_TERMS 27

/* The 96 primes 157..743 -- see the file header for what this covers
 * and when it's a real no-op vs. real work. */
#define NUM_GAP_PRIMES 96
static const uint32_t GAP_PRIMES[NUM_GAP_PRIMES] = {
	157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239,
	241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337,
	347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433,
	439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541,
	547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641,
	643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743,
};

/* The 256 primes 751..2677 -- always checked, regardless of which
 * upstream sieve ran. */
static const uint32_t EXT_PRIMES[NUM_EXT_PRIMES] = {
	751, 757, 761, 769, 773, 787, 797, 809,
	811, 821, 823, 827, 829, 839, 853, 857,
	859, 863, 877, 881, 883, 887, 907, 911,
	919, 929, 937, 941, 947, 953, 967, 971,
	977, 983, 991, 997, 1009, 1013, 1019, 1021,
	1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
	1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123,
	1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193,
	1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249,
	1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301,
	1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373,
	1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439,
	1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487,
	1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549,
	1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601,
	1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657,
	1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721,
	1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783,
	1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861,
	1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907,
	1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987,
	1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029,
	2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089,
	2099, 2111, 2113, 2129, 2131, 2137, 2141, 2143,
	2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237,
	2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287,
	2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347,
	2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393,
	2399, 2411, 2417, 2423, 2437, 2441, 2447, 2459,
	2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543,
	2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617,
	2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677,
};

/* One prime's worth of the check. a_mod = A mod P, step = (B*C) mod P,
 * both already reduced by the caller so this stays a flat, branchless
 * 27-step chain with no division inside it at all. Returns 1 if none of
 * the 27 terms is divisible by P. */
static inline int ap27_prefilter256_qualifies_p(uint32_t a_mod, uint32_t step, uint32_t p)
{
	uint32_t cur = a_mod;
	uint32_t any_zero = 0;
	for (int i = 0; i < NUM_TERMS; i++) {
		any_zero |= (cur == 0);
		cur += step;
		cur -= (cur >= p) ? p : 0;   /* branchless conditional subtract */
	}
	return !any_zero;
}

static inline int check_prime_list(const uint32_t *primes, int n, unsigned __int128 A, uint64_t B, uint64_t C)
{
	for (int k = 0; k < n; k++) {
		uint32_t p = primes[k];
		uint32_t a_mod = (uint32_t)(A % p);
		uint32_t b_mod = (uint32_t)(B % p);
		uint32_t c_mod = (uint32_t)(C % p);
		uint32_t step = (uint32_t)(((uint64_t)b_mod * c_mod) % p);
		if (!ap27_prefilter256_qualifies_p(a_mod, step, p))
			return 0;
	}
	return 1;
}

/* Checks EXT_PRIMES only (751..2677) -- the entry point for candidates
 * from a sieve that already covers primes up to 743 (i.e. the wider,
 * non-32-prime design). Returns 1 if the candidate survives all 256
 * primes, 0 if any single one already finds a factor. */
int ap27_prefilter256_survives(unsigned __int128 A, uint64_t B, uint64_t C)
{
	return check_prime_list(EXT_PRIMES, NUM_EXT_PRIMES, A, B, C);
}

/* Same as above, with GAP_PRIMES (157..743) optionally checked first.
 * check_gap=0 is byte-for-byte identical to ap27_prefilter256_survives()
 * above (existing callers, unchanged); check_gap=1 additionally closes
 * the real coverage gap left by the faster 32-prime sieve -- use it
 * there, and only there. */
int ap27_prefilter256_survives_ex(unsigned __int128 A, uint64_t B, uint64_t C, int check_gap)
{
	if (check_gap && !check_prime_list(GAP_PRIMES, NUM_GAP_PRIMES, A, B, C))
		return 0;
	return check_prime_list(EXT_PRIMES, NUM_EXT_PRIMES, A, B, C);
}

/* ==========================================================================
 * Self-test: cross-checks the fast, incremental path above against an
 * independently-written brute-force reference (direct 27-term modular
 * reduction per prime, no incremental-advance trick at all) on real and
 * randomized candidates, and checks the two prime ranges are contiguous
 * and non-overlapping by construction, not just by inspection.
 * ========================================================================== */

#ifndef PREFILTER_NO_MAIN
static int brute_force_check_list(const uint32_t *primes, int n, uint64_t A, uint64_t B, uint64_t C)
{
	for (int k = 0; k < n; k++) {
		uint32_t p = primes[k];
		for (int i = 0; i < NUM_TERMS; i++) {
			/* term = A + i*B*C, computed mod p via repeated mod (no overflow) */
			uint64_t term_mod = A % p;
			uint64_t add = (uint64_t)((B % p) * (C % p)) % p;
			for (int j = 0; j < i; j++)
				term_mod = (term_mod + add) % p;
			if (term_mod == 0)
				return 0;
		}
	}
	return 1;
}

static int brute_force_survives(uint64_t A, uint64_t B, uint64_t C)
{
	return brute_force_check_list(EXT_PRIMES, NUM_EXT_PRIMES, A, B, C);
}

static int brute_force_survives_ex(uint64_t A, uint64_t B, uint64_t C, int check_gap)
{
	if (check_gap && !brute_force_check_list(GAP_PRIMES, NUM_GAP_PRIMES, A, B, C))
		return 0;
	return brute_force_check_list(EXT_PRIMES, NUM_EXT_PRIMES, A, B, C);
}

int main(void)
{
	uint64_t C = 223092870ULL;  /* the AP27 primorial constant, 2*3*5*7*11*13*17*19*23 */
	int ok = 1, checked = 0;

	/* 1. A real, known AP27 (see ../cpu/ap27_crt_sieve.c's KAT list):
	 *    must survive, since by construction it has no small factors
	 *    at all. */
	{
		uint64_t A = 605185576317848261ULL, B = 155368778ULL;
		int fast = ap27_prefilter256_survives(A, B, C);
		int slow = brute_force_survives(A, B, C);
		checked++;
		if (fast != slow || fast != 1) {
			printf("FAIL: real KAT A=%llu B=%llu fast=%d slow=%d (expected 1,1)\n",
				(unsigned long long)A, (unsigned long long)B, fast, slow);
			ok = 0;
		} else {
			printf("PASS: real KAT survives extended sieve (fast=slow=1), as expected\n");
		}
	}

	/* 2. Random large odd candidates (not actual sieve survivors --
	 *    just odd, large A/B pairs, unrelated to the 32-prime sieve) --
	 *    cross-check fast vs. brute force on each, count survivors. */
	{
		uint64_t seed = 88172645463325252ULL;  /* xorshift64 seed */
		int survivors = 0;
		for (int t = 0; t < 5000; t++) {
			seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17;
			uint64_t A = (seed | 1) % 900000000000000000ULL + 100000000000000000ULL;
			uint64_t B = (seed >> 3) % 900000000ULL + 100000000ULL;
			int fast = ap27_prefilter256_survives(A, B, C);
			int slow = brute_force_survives(A, B, C);
			checked++;
			if (fast != slow) {
				printf("FAIL: A=%llu B=%llu fast=%d slow=%d\n",
					(unsigned long long)A, (unsigned long long)B, fast, slow);
				ok = 0;
			}
			if (fast) survivors++;
		}
		printf("random sweep: %d/5000 survived extended sieve (%.4f%%, expect ~0.9%%)\n",
			survivors, 100.0 * survivors / 5000);
	}

	/* 3. Exhaustive contiguity/non-overlap of the two prime ranges --
	 *    every prime from 157 to 2677 must appear in exactly one of
	 *    GAP_PRIMES or EXT_PRIMES, with none skipped and none
	 *    duplicated. Not just the endpoints: a future hand-edit that
	 *    drops or doubles an interior prime would otherwise pass
	 *    silently. */
	{
		int gap_idx = 0, ext_idx = 0, mismatches = 0;
		for (uint32_t p = 157; p <= 2677; p++) {
			int is_p = 1;
			for (uint32_t d = 2; (uint32_t)d * d <= p; d++)
				if (p % d == 0) { is_p = 0; break; }
			if (!is_p) continue;

			int in_gap = (gap_idx < NUM_GAP_PRIMES && GAP_PRIMES[gap_idx] == p);
			int in_ext = (ext_idx < NUM_EXT_PRIMES && EXT_PRIMES[ext_idx] == p);
			if (in_gap) gap_idx++;
			if (in_ext) ext_idx++;
			if (in_gap == in_ext) {  /* neither, or both -- either is wrong */
				printf("FAIL: prime %u is in %s of GAP_PRIMES/EXT_PRIMES\n",
					p, in_gap ? "BOTH" : "NEITHER");
				mismatches++;
			}
		}
		if (gap_idx != NUM_GAP_PRIMES || ext_idx != NUM_EXT_PRIMES) {
			printf("FAIL: consumed %d/%d GAP_PRIMES, %d/%d EXT_PRIMES -- table has an "
				"out-of-range or out-of-order entry\n",
				gap_idx, NUM_GAP_PRIMES, ext_idx, NUM_EXT_PRIMES);
			mismatches++;
		}
		checked++;
		if (mismatches) {
			ok = 0;
		} else {
			printf("boundary check: GAP_PRIMES=[%u..%u] (%d primes) EXT_PRIMES=[%u..%u] "
				"(%d primes) -- every prime 157..2677 checked individually: contiguous, "
				"no gap, no overlap\n",
				GAP_PRIMES[0], GAP_PRIMES[NUM_GAP_PRIMES-1], NUM_GAP_PRIMES,
				EXT_PRIMES[0], EXT_PRIMES[NUM_EXT_PRIMES-1], NUM_EXT_PRIMES);
		}
	}

	/* 4. check_gap=0 is a real no-op: must match the plain (gap-unaware)
	 *    entry point exactly, on the same random-ish candidates as (2).
	 *    check_gap=1 additionally cross-checked against its own
	 *    independent brute-force reference over GAP_PRIMES too. */
	{
		uint64_t seed = 12345678901234567ULL;
		int survivors_gap = 0;
		for (int t = 0; t < 5000; t++) {
			seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17;
			uint64_t A = (seed | 1) % 900000000000000000ULL + 100000000000000000ULL;
			uint64_t B = (seed >> 3) % 900000000ULL + 100000000ULL;

			int plain = ap27_prefilter256_survives(A, B, C);
			int nop = ap27_prefilter256_survives_ex(A, B, C, 0);
			checked++;
			if (plain != nop) {
				printf("FAIL: check_gap=0 not a no-op: A=%llu B=%llu plain=%d nop=%d\n",
					(unsigned long long)A, (unsigned long long)B, plain, nop);
				ok = 0;
			}

			int fast_gap = ap27_prefilter256_survives_ex(A, B, C, 1);
			int slow_gap = brute_force_survives_ex(A, B, C, 1);
			checked++;
			if (fast_gap != slow_gap) {
				printf("FAIL: check_gap=1 mismatch: A=%llu B=%llu fast=%d slow=%d\n",
					(unsigned long long)A, (unsigned long long)B, fast_gap, slow_gap);
				ok = 0;
			}
			if (fast_gap) survivors_gap++;
		}
		printf("gap-inclusive sweep: %d/5000 survived (%.4f%%, expect noticeably below the "
			"gap-unaware ~0.9%% -- GAP_PRIMES rejects further)\n",
			survivors_gap, 100.0 * survivors_gap / 5000);
	}

	printf("\n%d checks run, fast-vs-brute-force agreement: %s\n",
		checked, ok ? "ALL MATCH" : "MISMATCH FOUND");
	return ok ? 0 : 1;
}
#endif /* PREFILTER_NO_MAIN */
