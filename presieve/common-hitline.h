/*
 * common-hitline.h -- shared "B A" hit-line grammar + 128-bit decimal
 * helpers, used by every tool in this package that reads a search
 * driver's hit stream or a downstream survivor stream of it
 * (prefilter256_cli.c, ap27_log_decimator.c).
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * static-inline, single-header, "reuse as-is" style: no build step of
 * its own, just #include it from any .c file that needs the parsing/
 * printing helpers below.
 *
 * A is carried as a full unsigned __int128 throughout: real AP27
 * candidates are not bounded to 64 bits, so any tool that truncates A
 * before use risks a genuinely wrong result on a larger-A_lo run, not
 * just a cosmetic overflow.
 */
#ifndef AP27_COMMON_HITLINE_H
#define AP27_COMMON_HITLINE_H

#include <stdint.h>
#include <stdio.h>
#include <string.h>

/* Recognizes all 3 full-B line shapes this project's tools have emitted:
 *   "HIT B=<b> A=<a> ..."   (raw GPU/CPU driver stdout)
 *   "B=<b> A=<a>"           (an older survivor-stream shape)
 *   "<b> <a>"               (prefilter256_cli's own passthrough shape)
 * abuf must be >= 64 bytes. Returns 1 and fills *b/abuf on match, 0
 * otherwise (caller decides whether to skip or passthrough unmatched
 * lines). `a` may be plain decimal or `x`-prefixed hex (see
 * ap27_parse_u128_auto below) -- this function only captures the raw
 * token, it doesn't care which.
 */
static inline int ap27_hitline_parse(const char *line, unsigned long long *b, char *abuf)
{
	if (sscanf(line, "HIT B=%llu A=%63s", b, abuf) == 2) return 1;
	if (sscanf(line, "B=%llu A=%63s", b, abuf) == 2) return 1;
	if (sscanf(line, "%llu %63s", b, abuf) == 2) return 1;
	return 0;
}

/* decimal-ASCII -> unsigned __int128. No overflow/format validation --
 * matches this project's existing trust-the-input convention for its
 * own tool-to-tool pipes (the untrusted boundary is the Fermat check
 * itself, not this parse step).
 */
static inline unsigned __int128 ap27_parse_u128(const char *s)
{
	unsigned __int128 v = 0;
	while (*s >= '0' && *s <= '9') {
		v = v * 10 + (unsigned)(*s - '0');
		s++;
	}
	return v;
}

/* decimal-ASCII print of an unsigned __int128 to an arbitrary stream */
static inline void ap27_print_u128(FILE *out, unsigned __int128 v)
{
	char buf[48];
	int i = 48;
	if (v == 0) {
		fputc('0', out);
		return;
	}
	buf[--i] = 0;
	while (v > 0) {
		buf[--i] = '0' + (int)(v % 10);
		v /= 10;
	}
	fputs(&buf[i], out);
}

/* --- Compact hit-line format: hex A + truncated-B-suffix.
 * Motivation: at billions-of-lines scale (multi-day GPU/CPU runs across
 * many concurrent streams), decimal-formatting a value up to 128 bits
 * costs real per-line CPU (repeated division) in the write path --
 * enough to stall a many-thread run whose output feeds a live
 * compressed stream. Hex formatting is pure shift/mask. Full B is
 * logged exactly ONCE per stream, at start (ap27_print_b_header()
 * below is one ready-made way to do that; the CPU driver instead
 * prints its own ad-hoc header line with the same effect); every
 * subsequent hit line then carries only a 5-digit decimal suffix of B
 * -- enough to eyeball-disambiguate
 * interleaved concurrent streams in a merged log, NOT enough to
 * reconstruct B on its own (that needs the stream's own header line).
 * A is always full-precision hex, "x"-prefixed -- no "0x", just "x".
 *
 * Compact line shape: "HIT B=..<5 digits> A=x<hex digits> status=<word>"
 * status is always spelled out in full (never packed/coded) so it stays
 * trivially grep-able for any status word a downstream primality step
 * (not included in this repo) might set. The format is fully
 * self-describing (the "x" marks hex A, the ".." marks a B-suffix) so
 * a mixed log stream -- old-format and new-format lines interleaved,
 * or several streams at different digit widths -- parses unambiguously
 * line by line with no format flag or version marker needed.
 */

/* Full-B header line, written once per stream. */
static inline void ap27_print_b_header(FILE *out, unsigned long long b)
{
	fprintf(out, "B=%llu\n", b);
}

/* 5-digit decimal suffix of B, zero-padded so the field is always
 * exactly 5 characters wide (e.g. B=...00001 prints "..00001", not
 * "..1"). 5 digits keeps eyeball-collision risk low even at much higher
 * concurrent-stream counts than this project has run so far.
 */
static inline void ap27_print_b_suffix(FILE *out, unsigned long long b)
{
	fprintf(out, "..%05llu", (unsigned long long)(b % 100000ULL));
}

/* hex-ASCII (no "0x", just a leading "x") print of an unsigned __int128 */
static inline void ap27_print_u128_hex(FILE *out, unsigned __int128 v)
{
	static const char digits[] = "0123456789abcdef";
	char buf[40];
	int i = 40;
	buf[--i] = 0;
	if (v == 0) {
		buf[--i] = '0';
	} else {
		while (v > 0) {
			buf[--i] = digits[(unsigned)(v & 0xf)];
			v >>= 4;
		}
	}
	fputc('x', out);
	fputs(&buf[i], out);
}

/* decimal-or-hex ASCII -> unsigned __int128, auto-detected by a leading
 * 'x'/'X' (matches ap27_print_u128_hex's own output). Lets every
 * consumer transparently read both old (all-decimal) archives and new
 * compact (hex-A) logs through one call, with no format flag needed --
 * prefer this over the plain ap27_parse_u128() above in any new code,
 * and it's a safe drop-in replacement in existing parsers too (any
 * decimal-only archive still parses identically, since the 'x' check
 * only ever diverts genuinely hex-prefixed input).
 */
static inline unsigned __int128 ap27_parse_u128_auto(const char *s)
{
	if (*s == 'x' || *s == 'X') {
		unsigned __int128 v = 0;
		s++;
		for (;;) {
			char c = *s;
			int d;
			if (c >= '0' && c <= '9') d = c - '0';
			else if (c >= 'a' && c <= 'f') d = c - 'a' + 10;
			else if (c >= 'A' && c <= 'F') d = c - 'A' + 10;
			else break;
			v = (v << 4) | (unsigned)d;
			s++;
		}
		return v;
	}
	return ap27_parse_u128(s);
}

/* Recognizes the compact hit-line shape:
 *   "HIT B=..<5 digits> A=x<hex> status=<word>"
 * *bsuffix gets the 5-digit suffix (0..99999); abuf gets the raw A
 * token (pass to ap27_parse_u128_auto for the actual value). Returns 1
 * on match, 0 otherwise -- so a caller can try this after (or instead
 * of) ap27_hitline_parse() to handle a log stream with both line shapes
 * mixed together. Does NOT recover full B -- that requires the stream's
 * own "B=<decimal>" header line, read once up front by the caller and
 * matched against this suffix if full reconstruction is ever needed.
 */
static inline int ap27_hitline_parse_compact(const char *line, unsigned *bsuffix, char *abuf)
{
	return sscanf(line, "HIT B=..%5u A=%63s", bsuffix, abuf) == 2;
}

/* Extracts an arbitrary "<key><digits>" field from anywhere in the
 * line (e.g. key="batch=" finds "batch=1234" regardless of what other
 * fields surround it) -- lets the hit-line format grow new optional
 * fields over time (this project's drivers already have real examples:
 * batch=/pos=, the diagnostic indices a hit was computed from, printed
 * but not otherwise needed by most consumers) without needing a new
 * parser per field or per consumer. Returns 1 and fills *out if found,
 * 0 (leaving *out untouched) otherwise -- always safe against older
 * archives/lines that don't carry the field at all.
 */
static inline int ap27_hitline_field_u64(const char *line, const char *key, unsigned long long *out)
{
	const char *p = strstr(line, key);
	if (!p) return 0;
	return sscanf(p + strlen(key), "%llu", out) == 1;
}

/* --- Suffix -> full-B resolution map ---
 * A multi-worker run (several concurrent B's interleaved into ONE
 * stream) means a compact HIT line's 5-digit B-suffix alone isn't
 * enough to evaluate anything B-dependent downstream (e.g. a second-
 * stage sieve check) -- a single global "the stream's B" (fine for a
 * single-B GPU stream) is not even the right model there. This tracks
 * a small suffix->full-B map, populated by feeding it any full decimal
 * B a caller happens to observe (typically from a worker's own startup
 * header line, via ap27_hitline_field_u64 with key="B="), then
 * resolved against any subsequent compact HIT line sharing that
 * suffix.
 *
 * `static` (not just the functions `static inline`) so each
 * translation unit that includes this header gets its own private map
 * instance, exactly the "one map per running consumer" semantics
 * wanted here.
 *
 * Small, linear-scan table -- at most a handful of B's are ever
 * concurrently active (bounded by real thread counts), so this never
 * needs to be fast, just correct. Round-robin eviction if it somehow
 * fills (should not happen at any realistic concurrency, but never
 * grows unbounded either way).
 */
#define AP27_BMAP_SIZE 256
static unsigned ap27_bmap_suffix[AP27_BMAP_SIZE];
static unsigned long long ap27_bmap_full[AP27_BMAP_SIZE];
static int ap27_bmap_used[AP27_BMAP_SIZE];
static int ap27_bmap_evict_next = 0;

static inline void ap27_bmap_learn(unsigned long long full_b)
{
	unsigned suffix = (unsigned)(full_b % 100000ULL);
	for (int i = 0; i < AP27_BMAP_SIZE; i++) {
		if (ap27_bmap_used[i] && ap27_bmap_suffix[i] == suffix) {
			ap27_bmap_full[i] = full_b;
			return;
		}
	}
	for (int i = 0; i < AP27_BMAP_SIZE; i++) {
		if (!ap27_bmap_used[i]) {
			ap27_bmap_used[i] = 1;
			ap27_bmap_suffix[i] = suffix;
			ap27_bmap_full[i] = full_b;
			return;
		}
	}
	int i = ap27_bmap_evict_next;
	ap27_bmap_evict_next = (ap27_bmap_evict_next + 1) % AP27_BMAP_SIZE;
	ap27_bmap_suffix[i] = suffix;
	ap27_bmap_full[i] = full_b;
}

static inline int ap27_bmap_lookup(unsigned suffix, unsigned long long *out_full)
{
	for (int i = 0; i < AP27_BMAP_SIZE; i++) {
		if (ap27_bmap_used[i] && ap27_bmap_suffix[i] == suffix) {
			*out_full = ap27_bmap_full[i];
			return 1;
		}
	}
	return 0;
}

#endif /* AP27_COMMON_HITLINE_H */
