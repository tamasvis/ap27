/* common-cliargs.h -- shared dd-style "KEY=value" CLI argument parser.
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * Every multi-argument tool in this repository takes its numeric
 * arguments as "KEY=value" tokens, in ANY order, against a per-tool
 * table of expected keys and sane value bounds -- not a getopt()-style
 * flag parser, and not bare positional arguments. A mismatched or
 * unknown key fails immediately and exactly, and a value outside its
 * sane bound fails too, neither of which a bare positional parser can
 * check at all. A logged/history invocation line is also
 * self-documenting for free, since every token already says what it is.
 *
 * Internally, this header calls strlen() exactly where a length is not
 * yet known (argv entries, an interned string's first sighting), then
 * treats everything as a known-length byte buffer from then on --
 * memchr()/memcmp(), not strchr()/strcmp()/strncmp(). Chris Wellons:
 * "When I see strchr, strcmp or strncmp I wonder why you don't know
 * the lengths of your strings." (https://nullprogram.com/blog/2023/02/11/,
 * under "Strings")
 *
 * TODO: every failure path here is a hand-rolled `fprintf(stderr,
 * "FAIL: ...")` at the leaf where the bad input was actually detected,
 * deliberately, not a propagated/decoded error code: for a "should
 * not happen" parse failure like these, complaining loudly right
 * where it's detected beats building real error-code chaining and
 * decoding machinery just to report the same thing one call frame
 * later. If this project ever grows a shared "report this message to
 * stderr and return/pass through this code" helper, gated by its own
 * compile-time macro (defined: actually prints; undefined: silently
 * discards the message and just returns the code, so every call site
 * stays identical either way and reporting can be compiled out
 * entirely without touching them), these call sites are a reasonable
 * place to switch over to it. Not done now since a handful of
 * `fprintf` call sites don't justify a new dependency by themselves.
 */
#ifndef AP27_COMMON_CLIARGS_H
#define AP27_COMMON_CLIARGS_H

#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef enum { AP27_KV_U64, AP27_KV_STR } ap27_kv_type;

typedef struct {
	const char *key;
	ap27_kv_type type;
	/* AP27_KV_U64 fields: */
	uint64_t *u64_out; /* holds the parsed value on success, or whatever
	                    * the caller preset it to if optional and unseen */
	uint64_t min, max; /* inclusive; use UINT64_MAX for "no real ceiling" */
	/* AP27_KV_STR fields: */
	const char **str_out; /* on success, points into the static pool
	                        * below -- valid for the rest of the process.
	                        * If optional (required == 0), the caller
	                        * MUST preset *str_out to a real, non-NULL
	                        * default before calling: ap27_echo_kv_args()
	                        * always dereferences it, seen or not. */
	size_t str_maxlen;     /* excluding the terminating NUL */
	/* Common to both: */
	int required; /* 0: caller must preset the *_out field to the default */
	int seen;     /* internal; caller does not set this */
} ap27_kv_spec;

/* Backing storage for every string VALUE any ap27_kv_spec with type ==
 * AP27_KV_STR ever parses, across the whole process: one static
 * buffer, each accepted value appended as its own NUL-terminated run,
 * back to back, with a running "bytes used so far" offset -- the
 * simplest possible pool: no individual frees, no dynamic growth, and
 * a loud failure (not silent truncation) if it ever fills. A static,
 * bump-offset buffer like this is also a well-trodden, citable pattern
 * on its own terms -- see Chris Wellons, "Arena allocator tips and
 * tricks" (https://nullprogram.com/blog/2023/09/27/) for the general
 * technique (that post is about bump/arena allocation generally, not
 * string interning specifically; the interning/dedup behavior here is
 * this header's own addition on top of a plain arena).
 *
 * Deliberately CRUDE, and that is enough: a real CLI invocation here
 * is a handful of tokens, each at most a few dozen bytes, so matching
 * a repeated key's new value against its first occurrence by plain
 * memcmp -- no hashing, no dynamic growth -- costs nothing measurable
 * at this N.
 *
 * TODO, almost certainly never needed: if some future caller of this
 * header ever parses CLI args with enough string-valued keys/values
 * that scanning this pool linearly becomes measurably slow, or that
 * AP27_KV_POOL_BYTES below is too small, KEEP this pool as the flat
 * backing storage and add a hash table as an INDEX into it (key ->
 * pool offset) -- not a wholesale replacement of the pool, just a
 * lookup structure on top of it. See Wellons' follow-up, "An
 * easy-to-implement, arena-friendly hash map"
 * (https://nullprogram.com/blog/2023/09/30/), for a lightweight option
 * if that day ever comes. ap27_parse_kv_args() and ap27_echo_kv_args()
 * below never need to change either way, since both only ever see the
 * `const char *` this returns, never the pool itself. Don't make this
 * swap speculatively "just in case"; this crude version is more than
 * sufficient for every caller so far. */
#define AP27_KV_POOL_BYTES 4096
static char ap27_kv_pool[AP27_KV_POOL_BYTES];
static size_t ap27_kv_pool_used = 0;

static inline const char *ap27_kv_pool_intern(const char *s, size_t len)
{
	if (ap27_kv_pool_used + len + 1 > sizeof(ap27_kv_pool))
		return NULL; /* pool exhausted -- caller treats this as a parse failure */
	char *entry = ap27_kv_pool + ap27_kv_pool_used;
	memcpy(entry, s, len);
	entry[len] = '\0';
	ap27_kv_pool_used += len + 1;
	return entry;
}

/* Parses argv[start..argc) as KEY=value tokens, in any order, against
 * specs[0..nspecs). Strict: an unrecognized key, a malformed/out-of-
 * range/too-long value, or a missing required key all fail. A key
 * repeated with the SAME value (byte-for-byte, for a string spec) is
 * accepted as a harmless no-op (so "A=5 A=5" or "MODE=fast MODE=fast"
 * is fine); repeated with a DIFFERENT value, it fails -- that is much
 * more likely to be a copy-paste mistake than an intentional override,
 * and this parser has no "last one wins" convention to make such a
 * mistake silent. Prints a usage message to stderr and returns 0 on
 * any failure. usage should already end in '\n'. */
static inline int ap27_parse_kv_args(int argc, char **argv, int start,
                                      ap27_kv_spec *specs, int nspecs,
                                      const char *usage)
{
	for (int s = 0; s < nspecs; s++) specs[s].seen = 0;
	for (int i = start; i < argc; i++) {
		size_t toklen = strlen(argv[i]);
		const char *eq = (const char *)memchr(argv[i], '=', toklen);
		if (!eq) {
			fprintf(stderr, "FAIL: expected KEY=value, got \"%s\"\n%s", argv[i], usage);
			return 0;
		}
		size_t klen = (size_t)(eq - argv[i]);
		int s;
		for (s = 0; s < nspecs; s++)
			if (strlen(specs[s].key) == klen && memcmp(argv[i], specs[s].key, klen) == 0)
				break;
		if (s == nspecs) {
			fprintf(stderr, "FAIL: unrecognized key in \"%s\"\n%s", argv[i], usage);
			return 0;
		}
		const char *valstr = eq + 1;

		if (specs[s].type == AP27_KV_STR) {
			size_t vlen = strlen(valstr);
			if (vlen == 0 || vlen > specs[s].str_maxlen) {
				fprintf(stderr, "FAIL: %s=%s is empty or longer than %zu characters\n%s",
					specs[s].key, valstr, specs[s].str_maxlen, usage);
				return 0;
			}
			if (specs[s].seen) {
				size_t prevlen = strlen(*specs[s].str_out);
				if (prevlen != vlen || memcmp(*specs[s].str_out, valstr, vlen) != 0) {
					fprintf(stderr, "FAIL: %s given twice with different values (\"%s\", then \"%s\")\n%s",
						specs[s].key, *specs[s].str_out, valstr, usage);
					return 0;
				}
			}
			if (!specs[s].seen) {
				const char *interned = ap27_kv_pool_intern(valstr, vlen);
				if (!interned) {
					fprintf(stderr, "FAIL: internal string pool exhausted parsing %s=...\n%s",
						specs[s].key, usage);
					return 0;
				}
				*specs[s].str_out = interned;
			}
			specs[s].seen = 1;
			continue;
		}

		/* strtoull() on its own accepts leading whitespace and a
		 * leading '+'/'-' (so "-1" parses as a huge positive value via
		 * wraparound, not an error), and silently clamps on overflow
		 * (returning ULLONG_MAX with errno=ERANGE, which must be
		 * checked explicitly -- it is not reflected in *end). Requiring
		 * the first character to already be a digit, and checking
		 * errno after the call, closes both: "FAIL: expected KEY=value"
		 * exists specifically so a malformed value can't be silently
		 * misread as a different, valid-looking one. */
		if (valstr[0] < '0' || valstr[0] > '9') {
			fprintf(stderr, "FAIL: %s=%s is not a plain decimal integer\n%s",
				specs[s].key, valstr, usage);
			return 0;
		}
		errno = 0;
		char *end;
		uint64_t v = strtoull(valstr, &end, 10);
		if (*end != '\0' || errno == ERANGE) {
			fprintf(stderr, "FAIL: %s=%s is not a plain decimal integer\n%s",
				specs[s].key, valstr, usage);
			return 0;
		}
		if (v < specs[s].min || v > specs[s].max) {
			fprintf(stderr, "FAIL: %s=%llu is outside the allowed range [%llu, %llu]\n%s",
				specs[s].key, (unsigned long long)v,
				(unsigned long long)specs[s].min, (unsigned long long)specs[s].max, usage);
			return 0;
		}
		if (specs[s].seen && *specs[s].u64_out != v) {
			fprintf(stderr, "FAIL: %s given twice with different values (%llu, then %llu)\n%s",
				specs[s].key, (unsigned long long)*specs[s].u64_out, (unsigned long long)v, usage);
			return 0;
		}
		*specs[s].u64_out = v;
		specs[s].seen = 1;
	}
	for (int s = 0; s < nspecs; s++) {
		if (specs[s].required && !specs[s].seen) {
			fprintf(stderr, "FAIL: missing required %s=...\n%s", specs[s].key, usage);
			return 0;
		}
	}
	return 1;
}

/* Echoes the parsed values as "prog: KEY=value KEY=value ..." to
 * stderr -- call right after a successful ap27_parse_kv_args() so the
 * invocation is self-documenting in any log/history that captures
 * stderr, independent of whatever each tool's own banner prints. */
static inline void ap27_echo_kv_args(const char *prog, ap27_kv_spec *specs, int nspecs)
{
	fprintf(stderr, "%s:", prog);
	for (int s = 0; s < nspecs; s++) {
		if (specs[s].type == AP27_KV_STR)
			fprintf(stderr, " %s=%s", specs[s].key, *specs[s].str_out);
		else
			fprintf(stderr, " %s=%llu", specs[s].key, (unsigned long long)*specs[s].u64_out);
	}
	fprintf(stderr, "\n");
}

#endif
