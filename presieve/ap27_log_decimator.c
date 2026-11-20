/*
 * ap27_log_decimator.c -- streaming second-stage filter for a search
 * driver's raw hit stream. Meant to sit between the driver and the
 * compressor: `search | ap27_log_decimator | zstd ...`.
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * A 32-prime sieve survivor is common enough (~3.48e-6 of all
 * candidates) that logging every single one becomes a real disk-
 * volume constraint on long unattended runs, well before it becomes
 * anything else. But only ~5.18e-6 of those survivors also pass the
 * second-stage sieve in this same directory (ap27_prefilter256.c,
 * GAP_PRIMES union EXT_PRIMES, 352 primes 157..2677) -- a ~193,060x
 * reduction on top of the 32-prime sieve's own filtering. There is no
 * reason to persist the full first-stage stream when a much smaller,
 * already-correct downstream filter already exists; this tool applies
 * it inline, before compression, rather than as an offline pass over
 * an already-huge archive.
 *
 * What gets logged, and why:
 *   - every second-stage survivor: always, in full. Rare, and by
 *     definition worth a closer look.
 *   - every new widest-gap-yet among first-stage survivors (top K):
 *     always, in full, regardless of that candidate's own second-
 *     stage outcome. Gap-worthiness is a property of the first-stage
 *     arrival process itself, independent of what the second-stage
 *     filter thinks of that specific candidate.
 *   - a small, configurable sample of ordinary first-stage hits: just
 *     enough to spot-check the observed rate against the theoretical
 *     one, without paying for every single one.
 *   - anything that isn't a HIT line (headers, progress, DONE) is
 *     passed through completely unchanged.
 *
 * Gap-tracking has to be robust to a real property of the parallel
 * GPU sieve: hits are not guaranteed to arrive in strictly increasing-
 * A order (concurrent invocations racing on the atomic hit-count
 * append can interleave nearby A's within a chunk). A running MAXIMUM
 * A, not the previous line's A, is what a gap is measured against --
 * a late-arriving, already-covered value is folded into the total
 * count and otherwise ignored, never misread as a large negative gap.
 *
 * Output uses PLAIN, full-precision decimal B/A: a compact hex-A/
 * B-suffix shape exists upstream specifically to survive a many-
 * billion-line-per-run firehose without paying per-line decimal-
 * formatting CPU or storage cost, but that justification doesn't apply
 * to THIS tool's own output once it has already cut that volume by
 * ~193,060x -- there's no remaining reason to also pay the compact
 * format's reduced-readability cost on what's left; full plain decimal
 * (matching testvectors/hitfilter_kat_survivor1.txt's format) is more
 * directly usable downstream (eyeballing, or feeding to a real
 * Fermat/Miller-Rabin confirmation step) for a comparatively tiny
 * remaining line count. Compression (zstd
 * downstream of this tool) is a SEPARATE, independent decision from
 * this -- it stays on regardless, for pipeline consistency, even
 * though its own marginal benefit shrinks once decimation has already
 * done most of the volume reduction.
 *
 * INPUT, by contrast, still has to cope with the compact shape: the
 * search drivers upstream of this tool (e.g. ap27_crt_sieve --timed)
 * emit compact HIT lines (only a 5-digit B-suffix per line) for their
 * OWN, still-real reason (undecimated per-line cost at the driver's own
 * write-path), and a real multi-worker run (several concurrent B's
 * interleaved into ONE stream) needs the FULL B to evaluate the
 * second-stage filter correctly. A small suffix->full-B
 * map handles this: populated by matching each worker's own startup
 * header line ("B=<full decimal>", found anywhere in a non-HIT line
 * via ap27_hitline_field_u64) as it passes through, then resolved
 * against any subsequent compact HIT line sharing that same suffix. A
 * suffix with no known full B yet (e.g. a truncated stream missing its
 * own header) passes through unchanged rather than guessing.
 *
 * Usage: <search driver> | ap27_log_decimator [SAMPLE_EVERY] [--verify] | zstd ...
 *
 * --verify: cross-check every incoming first-stage HIT line's A against
 * all primes up to 151 (the full set the upstream 32-prime sieve is
 * supposed to enforce, via its own tables/gates) and abort loudly,
 * nonzero exit, the moment one fails -- a failure here means the
 * upstream sieve is not rejecting what it claims to, which no existing
 * check in this project's `make check` can otherwise detect (a weakened
 * or disabled per-prime table still passes every KAT, since KATs only
 * test for a known true positive, never that everything else is
 * correctly rejected). Off by default: the check re-verifies primes the
 * fast ROM-table sieve already specialized tables for, so it costs real,
 * avoidable CPU on a trusted production stream -- opt in only when
 * actually auditing a build or a suspicious run. */

#define PREFILTER_NO_MAIN
#include "ap27_prefilter256.c"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common-hitline.h"
#include "../common-cliargs.h"

#define TOPK 20
#define DEFAULT_SAMPLE_EVERY 100000
#define SUMMARY_EVERY 1000000
/* One named size for every "a real hit/progress line always fits in
 * this" buffer, rather than the same-looking magic number picked
 * separately in several places. A real compact hit line (B suffix +
 * up to 128-bit hex A + status + batch=/pos=) fits well under 128
 * bytes; this leaves real headroom. */
#define LINE_BIG_ENOUGH 256

static const uint64_t C_VAL_LOCAL = 223092870ull;

/* Every prime up to 151, for --verify. The upstream 32-prime sieve
 * enforces this exact set, just split across three mechanisms (2,3 via
 * a mod-6 gate; 29,31 via the mod-5394 gate; the rest via per-prime ROM
 * tables, 5..23 CRT-merged) -- check_prime_list() below doesn't care how
 * a given prime is normally enforced, only whether it actually is. */
#define NUM_VERIFY_PRIMES 36
static const uint32_t VERIFY_PRIMES[NUM_VERIFY_PRIMES] = {
	2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
	67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137,
	139, 149, 151
};

typedef struct {
	unsigned __int128 gap;
	char prev_line[LINE_BIG_ENOUGH];
	char cur_line[LINE_BIG_ENOUGH];
} gap_entry;

static gap_entry topk[TOPK];
static int topk_n = 0;

/* Per-stream running-max-A table, keyed by the driver's own `id=` field:
 * distinct concurrent workers can share the same B while scanning
 * different A ranges, so B alone isn't a safe gap-tracking key. Lines
 * with no `id=` field (the common case here -- the GPU
 * driver doesn't currently emit one, since it only ever runs a single
 * stream per process) all share stream id 0, which is correct as long
 * as that one id is only ever used by a single real stream at a time.
 *
 * Small, linear-scan table -- at most a handful of concurrent streams
 * ever exist, so this never needs to be fast, just correct. Round-
 * robin eviction if it somehow fills. */
#define STREAM_TABLE_ELEMS 64
typedef struct {
	unsigned id;
	unsigned __int128 running_max_A;
	int have_max;
	int used;
	char last_raw_line[LINE_BIG_ENOUGH];
} stream_state;
static stream_state stream_table[STREAM_TABLE_ELEMS];
static int stream_evict_next = 0;

static stream_state *stream_state_get(unsigned id)
{
	for (int i = 0; i < STREAM_TABLE_ELEMS; i++)
		if (stream_table[i].used && stream_table[i].id == id)
			return &stream_table[i];
	for (int i = 0; i < STREAM_TABLE_ELEMS; i++) {
		if (!stream_table[i].used) {
			stream_table[i].used = 1;
			stream_table[i].id = id;
			stream_table[i].have_max = 0;
			stream_table[i].running_max_A = 0;
			stream_table[i].last_raw_line[0] = '\0';
			return &stream_table[i];
		}
	}
	int i = stream_evict_next;
	stream_evict_next = (stream_evict_next + 1) % STREAM_TABLE_ELEMS;
	stream_table[i].used = 1;
	stream_table[i].id = id;
	stream_table[i].have_max = 0;
	stream_table[i].running_max_A = 0;
	stream_table[i].last_raw_line[0] = '\0';
	return &stream_table[i];
}

/* Inserts (gap, prev, cur) into the sorted (descending) top-K list if
 * it qualifies. Returns 1 if it was inserted (a new record among the
 * K widest seen so far), 0 otherwise. */
static int topk_maybe_insert(unsigned __int128 gap, const char *prev, const char *cur)
{
	if (topk_n == TOPK && gap <= topk[TOPK - 1].gap)
		return 0;
	int pos = (topk_n < TOPK) ? topk_n : TOPK - 1;
	if (topk_n < TOPK) topk_n++;
	while (pos > 0 && topk[pos - 1].gap < gap) {
		topk[pos] = topk[pos - 1];
		pos--;
	}
	topk[pos].gap = gap;
	snprintf(topk[pos].prev_line, sizeof(topk[pos].prev_line), "%s", prev);
	snprintf(topk[pos].cur_line, sizeof(topk[pos].cur_line), "%s", cur);
	return 1;
}

/* Plain, full-precision output -- see this file's own header comment
 * for why: the compact shape's justification doesn't survive
 * decimation, and readability/downstream-tool-friendliness wins once
 * the line count is this small. */
static void print_plain(unsigned long long B, unsigned __int128 A, const char *tag)
{
	printf("HIT B=%llu A=", B);
	ap27_print_u128(stdout, A);
	printf(" status=unchecked stage2=%s\n", tag);
}

/* Same shape as print_plain(), plus an extra additive-safe `gap=` field
 * (same spirit as stage2=/batch=/pos=) carrying the EXACT gap value at
 * the moment it was computed -- this is the only place that value is
 * ever known; the end-of-run "FINAL: top K widest gaps" stderr report
 * is a convenience summary, not the source of truth, and a caller
 * redirecting this tool's stderr away (e.g. to /dev/null) loses it
 * entirely. Putting the value directly on the stdout line that's
 * already being archived means the true top-K is always exactly
 * reconstructable later by a plain sort over `gap=`, independent of
 * what happens to stderr or whether this process's own FINAL summary
 * ever reaches anywhere durable. */
static void print_gap_record(unsigned long long B, unsigned __int128 A, unsigned __int128 gap)
{
	printf("HIT B=%llu A=", B);
	ap27_print_u128(stdout, A);
	printf(" status=unchecked stage2=gap_record gap=");
	ap27_print_u128(stdout, gap);
	printf("\n");
}

int main(int argc, char **argv)
{
	/* --verify is a bare flag (no value), so it's pulled out of argv
	 * before the rest goes through the shared KEY=value parser (which
	 * -- same reasoning as the sieve drivers' B=/A=/SECONDS= args --
	 * rejects a typo like --verfy outright instead of silently parsing
	 * it as 0 and leaving verify_mode off with no warning, a real trap
	 * for an audit flag). argc is always tiny here, so a plain
	 * fixed-size local array for the filtered-down KEY=value tokens
	 * needs no bound checking beyond argc itself. */
	int verify_mode = 0;
	char *kvargv[argc > 0 ? argc : 1];
	int nkv = 0;
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "--verify") == 0)
			verify_mode = 1;
		else
			kvargv[nkv++] = argv[i];
	}
	char usagebuf[128];
	snprintf(usagebuf, sizeof(usagebuf), "usage: %s [--verify] [SAMPLE_EVERY=<n>]\n", argv[0]);
	uint64_t sample_every_u64 = (uint64_t)DEFAULT_SAMPLE_EVERY;
	ap27_kv_spec specs[] = {
		{ .key = "SAMPLE_EVERY", .type = AP27_KV_U64, .u64_out = &sample_every_u64,
		  .min = 1, .max = 1000000000ull },
	};
	if (!ap27_parse_kv_args(nkv, kvargv, 0, specs, 1, usagebuf))
		return 1;
	ap27_echo_kv_args("ap27_log_decimator", specs, 1);
	fprintf(stderr, "ap27_log_decimator: --verify=%s\n", verify_mode ? "on" : "off");
	long sample_every = (long)sample_every_u64;

	/* Config banner on STDOUT, not stderr -- deliberately, so it's part
	 * of the same archived stream as the HIT lines themselves and
	 * survives regardless of how a caller handles this process's
	 * stderr. Lets a reader decode a bare log file's sampling
	 * ratio/TOPK without needing outside context or memory of what
	 * this run was launched with. */
	printf("[log_decimator] config: sample_every=%ld topk=%d verify=%s\n",
		sample_every, TOPK, verify_mode ? "on" : "off");
	fflush(stdout);

	unsigned long long B;
	char abuf[64];
	char line[LINE_BIG_ENOUGH];
	long total = 0, survived = 0, sampled = 0, records = 0;

	while (fgets(line, sizeof(line), stdin)) {
		unsigned bsuffix;
		char abuf_compact[64];
		if (ap27_hitline_parse_compact(line, &bsuffix, abuf_compact)) {
			/* multi-worker CPU stream shape: resolve the suffix against
			 * whichever worker's header line established it. Unknown
			 * suffix (e.g. a truncated stream missing its own header)
			 * passes through unchanged rather than guessing which B. */
			if (!ap27_bmap_lookup(bsuffix, &B)) {
				fputs(line, stdout);
				continue;
			}
			memcpy(abuf, abuf_compact, sizeof(abuf));
		} else if (!ap27_hitline_parse(line, &B, abuf)) {
			unsigned long long learned;
			if (ap27_hitline_field_u64(line, "B=", &learned))
				ap27_bmap_learn(learned);
			fputs(line, stdout);
			continue;
		}
		total++;
		unsigned __int128 A = ap27_parse_u128_auto(abuf);

		if (verify_mode && !check_prime_list(VERIFY_PRIMES, NUM_VERIFY_PRIMES, A, B, C_VAL_LOCAL)) {
			fflush(stdout);
			fprintf(stderr,
				"[log_decimator] FAIL: --verify: line %ld claims to be a "
				"32-prime-sieve survivor but is divisible by a term <= 151 "
				"(B=%llu A=", total, B);
			ap27_print_u128(stderr, A);
			fprintf(stderr, ") -- the upstream sieve is not rejecting "
				"everything it claims to; offending line:\n%s", line);
			return 1;
		}

		/* gap-tracking, per-stream -- see stream_state_get()'s own
		 * comment for why B alone isn't a safe key. Lines with no id=
		 * field (the common case for this tool -- see that comment)
		 * all share stream id 0. */
		unsigned long long id64 = 0;
		ap27_hitline_field_u64(line, "id=", &id64);
		stream_state *st = stream_state_get((unsigned)id64);

		if (st->have_max && A > st->running_max_A) {
			unsigned __int128 gap = A - st->running_max_A;
			if (topk_maybe_insert(gap, st->last_raw_line, line)) {
				records++;
				print_gap_record(B, A, gap);
			}
		}
		if (!st->have_max || A > st->running_max_A) {
			st->running_max_A = A;
			st->have_max = 1;
		}
		snprintf(st->last_raw_line, sizeof(st->last_raw_line), "%s", line);

		if (ap27_prefilter256_survives_ex(A, B, C_VAL_LOCAL, 1)) {
			survived++;
			print_plain(B, A, "survivor");
		} else if (total % sample_every == 0) {
			sampled++;
			print_plain(B, A, "sampled");
		}

		if (total % SUMMARY_EVERY == 0) {
			fprintf(stderr,
				"[log_decimator] total=%ld survived=%ld (%.6f%%, expect ~0.000518%%) "
				"sampled=%ld gap_records=%ld\n",
				total, survived, total ? 100.0 * survived / total : 0.0,
				sampled, records);
			fflush(stderr);
		}
	}

	fprintf(stderr, "[log_decimator] FINAL: total=%ld survived=%ld (%.6f%%) sampled=%ld gap_records=%ld\n",
		total, survived, total ? 100.0 * survived / total : 0.0, sampled, records);
	fprintf(stderr, "[log_decimator] top %d widest gaps:\n", topk_n);
	for (int i = 0; i < topk_n; i++) {
		fprintf(stderr, "  #%d gap=", i + 1);
		ap27_print_u128(stderr, topk[i].gap);
		fprintf(stderr, "\n    prev: %s    cur:  %s", topk[i].prev_line, topk[i].cur_line);
	}
	if (verify_mode && total == 0) {
		fprintf(stderr, "[log_decimator] FAIL: --verify audited 0 HIT lines -- "
			"nothing was checked (did the upstream sieve fail or print nothing?)\n");
		return 1;
	}
	return 0;
}
