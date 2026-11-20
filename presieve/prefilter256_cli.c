/* prefilter256_cli.c -- streaming filter, reads sieve survivors on stdin,
 * writes presieve survivors (candidates ready for a Fermat/Miller-Rabin
 * test) to stdout.
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * Accepts every hit-line shape this project's drivers emit: the GPU
 * drivers' full-decimal "B=<b> A=<a>"/"HIT B=<b> A=<a> ..." lines, the
 * CPU driver's compact "HIT B=..<suffix> A=x<hex> ..." lines (resolved
 * against a suffix->full-B map learned from each stream's own header
 * line, same as ap27_log_decimator.c), or plain "<b> <a>"; anything
 * else is passed to stderr, not emitted as a candidate. Emits surviving
 * "<b> <a>" pairs, one per line, with A always decimal regardless of
 * which input shape it was read from (a downstream Fermat/MR stage
 * gets one consistent format, not hex for one stream and decimal for
 * another).
 *
 * GAP_PRIMES (157..743) are checked by default: this file's companion
 * drivers (../cpu/, ../gpu/) only sieve up to 151, so that range is a
 * real, otherwise-unchecked gap for their output. Pass --no-gap for
 * input from a wider sieve that already covers primes up to 743 (where
 * checking again would be a real no-op, just wasted work).
 *
 * Build: make (from the repo root) -- see the Makefile for the exact
 *        flags (deliberately -O3 -march=native, not just the project
 *        default -O2, for this specific tool), rather than duplicating
 *        them here where they could drift out of sync.
 * Run:   ./prefilter256_cli [--no-gap] < hits.txt > survivors.txt
 */
#define PREFILTER_NO_MAIN
#include "ap27_prefilter256.c"
#include <stdio.h>
#include <string.h>
#include "common-hitline.h"

int main(int argc, char **argv)
{
	int check_gap = !(argc > 1 && strcmp(argv[1], "--no-gap") == 0);
	uint64_t C = 223092870ULL;
	unsigned long long B;
	char abuf[64];
	char line[256];
	long total = 0, survived = 0;
	while (fgets(line, sizeof(line), stdin)) {
		unsigned bsuffix;
		char abuf_compact[64];
		if (ap27_hitline_parse_compact(line, &bsuffix, abuf_compact)) {
			if (!ap27_bmap_lookup(bsuffix, &B)) {
				fputs(line, stderr);
				continue;
			}
			memcpy(abuf, abuf_compact, sizeof(abuf));
		} else if (!ap27_hitline_parse(line, &B, abuf)) {
			unsigned long long learned;
			if (ap27_hitline_field_u64(line, "B=", &learned))
				ap27_bmap_learn(learned);
			fputs(line, stderr);
			continue;
		}
		total++;
		unsigned __int128 A = ap27_parse_u128_auto(abuf);
		if (ap27_prefilter256_survives_ex(A, B, C, check_gap)) {
			printf("%llu ", B);
			ap27_print_u128(stdout, A);
			putchar('\n');
			survived++;
		}
	}
	fprintf(stderr, "prefilter256%s: %ld/%ld survived (%.4f%%)\n",
		check_gap ? "+gap" : "", survived, total, total ? 100.0*survived/total : 0.0);
	return 0;
}
