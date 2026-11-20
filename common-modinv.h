/*
 * common-modinv.h -- shared modular-inverse helper for the AP27 CRT/
 * ROM-batch sieve's per-prime table construction.
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * `cpu/` and `gpu/` are sibling directories of this same package and
 * both need this exact function for their ROM-table construction, so
 * it lives once at the package root and both sides reach it via
 * `#include "../common-modinv.h"` -- true single-source-of-truth,
 * rather than two copies that need to be kept in sync by hand.
 * `static inline` so this can be #included directly into multiple .c
 * files without violating the one-definition rule.
 *
 * Usage: #include "../common-modinv.h" from cpu/ or gpu/, call
 * modinv_u64(a, m) for any prime modulus m. Pulls in <stdint.h>,
 * <stdio.h>, <stdlib.h> itself -- a caller doesn't need to include
 * them first.
 */

#ifndef AP27_COMMON_MODINV_H
#define AP27_COMMON_MODINV_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/* Extended-Euclid modular inverse of a mod m. Real inverse only exists
 * when gcd(a,m)==1 -- for the prime moduli this project always calls
 * this with, that means whenever a%m==0 (the only way gcd(a,m)>1 for
 * prime m). The plain extended-Euclid loop below doesn't detect that on
 * its own and can return a garbage-but-plausible value instead of
 * erroring; it is intended to be simple, not fast. m=0 input
 * falls through to t=0 as an (incorrect) answer.
 *
 * Post-condition check: verify (a*result) mod m == 1 before ever
 * returning it, rather than trusting the loop's output. Callers are
 * welcome to add their own additional checks on top.
 */
static inline uint64_t modinv_u64(uint64_t a, uint64_t m)
{
	int64_t t = 0, newt = 1;
	int64_t r = (int64_t)m, newr = (int64_t)a;

	while (newr) {
		int64_t q = r / newr;
		int64_t tmp_t = t - q * newt; t = newt; newt = tmp_t;
		int64_t tmp_r = r - q * newr; r = newr; newr = tmp_r;
	}
	if (t < 0) t += (int64_t)m;

	uint64_t check_product = (uint64_t)(((__int128)a * (__int128)t) % (__int128)m);
	uint64_t expect_one = (uint64_t)(1 % m);
	if (check_product != expect_one) {
		fprintf(stderr, "FAIL: modinv_u64(%llu, %llu) = %llu is not a real "
			"inverse (a*result mod m = %llu, expected %llu) -- gcd(a,m)>1, no "
			"inverse exists; this should have been caught by a caller-side "
			"B-rejection guard before ever reaching here\n",
			(unsigned long long)a, (unsigned long long)m, (unsigned long long)t,
			(unsigned long long)check_product, (unsigned long long)expect_one);
		exit(1);
	}

	return (uint64_t)t;
}

#endif /* AP27_COMMON_MODINV_H */
