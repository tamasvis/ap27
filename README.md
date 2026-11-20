# AP27 CRT/ROM-batch prime-progression sieve

Visegrady, Tamas ```<tamas@visegrady.ch>```

*Assisted by Claude (claude-sonnet-5-5).*

We describe a sieve for AP27 search---27 primes in arithmetic
progression---built around one idea: a candidate's sieve state is its
residue tuple modulo a fixed set of primes---a CRT (Chinese Remainder
Theorem) representation---and the per-prime "is this candidate provably
composite?" answer for the next 512 candidates at once is a function of
*a precomputed bitmask, indexed purely by residue*, never by position.

Since we only check fixed-sized, constant bitmasks ANDed together,
never divide at runtime, and need just a few dozen bytes of updated
state, our design's hot path does less work per candidate than one that
recalculates wider bitmasks as the search advances. The small state also
simplifies parallelizing the search even on smaller machines.

## The problem

An AP27 is 27 primes in arithmetic progression: `A, A+D, A+2D,
..., A+26D`, where `D = B*C` for a fixed constant `C = 223092870`
(`2*3*5*7*11*13*17*19*23`, the primorial of the primes up to 23) and a
free parameter `B`. Searching efficiently means sweeping `A` for a fixed
`B` and rejecting many candidates cheaply, without ever computing the
27 actual terms---a candidate is not an AP27 if any term provably has a
small prime as a factor.

## The design

Most integers are not worth checking: potential candidate positions
follow a fixed period of 5394 (`= 6*29*31`), determined once per `B`;
there are always 16 possible offsets[^1] in each period. The next 512
candidate offsets are fixed in 32 adjacent mod-5394 periods. For any
prime `P`, "does offset `s` (one of the next 512) indicate that the
`(A+s,B)` AP27-sequence has a term divisible by `P`?" has a fixed
512-bit answer vector which depends **only on `A mod P`** and **fixed
tables specialized from `B`** -- these are the "ROM" in this design's
name: precomputed, read-only lookup tables, built once per `B` and
not modified mid-search. (In other words, we only need a CRT-based
representation of `A,` not `A` itself.)

We scan `A` candidates for a fixed `B` and precompute a 512-bit
eligibility bitmask for every possible `mod P` residue---fixed for
any specific `B`, specialized during initialization (below, `C#23` is
the reference script's own shorthand for `C`, the primorial of the
primes up to 23):

```
B=1234567

The 16 possible mod-6*29*31 offsets (function of B) are:
[11, 107, 397, 707, 803, 1499, ..., 4709, 5095]
512-bit bitmasks store offsets for the next 32 full 6*29*31 periods.

'A' candidates are advanced to the first A with A mod 5394 == 11, then
we use offsets relative to that start:
[0, 96, 386, 696, 792, 1488, ...

For each small prime 'P', we evaluate whether any A + B * C#23 * i (i in
[0..26]) is a multiple of P at each of the 512 next tabulated offsets
(here showing the start with current position at A = 1 mod P; cand%P
is that specific candidate's own residue, not the batch start's):

P=53    (B=38 mod P; C#23 = 23 mod P; B * C#23 = 26 mod P)
BIT[  0][     +0] cand%P=1   ## A=53*k+1 -> 53|A + (B=1234567) * C#23 * 2
BIT[  1][    +96] cand%P=44
BIT[  2][   +386] cand%P=16
BIT[  3][   +696] cand%P=8   ## A=53*k+8 -> 53|A + (B=1234567) * C#23 * 16
BIT[  4][   +792] cand%P=51
BIT[  5][  +1488] cand%P=5   ## A=53*k+5 -> 53|A + (B=1234567) * C#23 * 10
BIT[  6][  +2184] cand%P=12  ## A=53*k+12 -> 53|A + (B=1234567) * C#23 * 24
BIT[  7][  +2204] cand%P=32  ## A=53*k+32 -> 53|A + (B=1234567) * C#23 * 11
...
BIT[504][+170114] cand%P=38  ## A=53*k+38 -> 53|A + (B=1234567) * C#23 * 23
BIT[505][+170810] cand%P=45
BIT[506][+170906] cand%P=35  ## A=53*k+35 -> 53|A + (B=1234567) * C#23 * 17
BIT[507][+171216] cand%P=27  ## A=53*k+27 -> 53|A + (B=1234567) * C#23
BIT[508][+171506] cand%P=52
BIT[509][+171602] cand%P=42
BIT[510][+171912] cand%P=34  ## A=53*k+34 -> 53|A + (B=1234567) * C#23 * 15
BIT[511][+172298] cand%P=49

The above 'possible' offsets (bit=1), when packed left-to-right as 512
bits, form (real output, not abbreviated):
## BITMASK[P=53][248/512]=0x6845953a62c59d1b62e68d1972a6885937a2e85d97b26845953a62c49d1b72e68d1973a6c85917a2e84d97b268c59d3a62c49d1972e68d5977a6c85917a2684d
(248 are possibly AP27's out of 512)
```
Reference code in ```ap27sample.py``` shows full construction
details.

Each batched-sieve step fetches a bitmask per prime `P`, ANDs them
all together; checks any offset for a possible AP27 where the ANDed
bit is `1`. Moving on to the next batch advances `A` with the fixed
increment of `32*6*29*31,` which requires only one modular addition per
P prime---update an index---using a precomputed table.

[^1]: modulo 29 and 31 restrictions define most possible offsets, which
are then filtered to skip even or multiple-of-3 positions (the `A mod 6`
gate for 2 and 3).

This design's basic "check many candidates via a bitmask"
framing was originally inspired by the reference PrimeGrid
AP26 client's own AVX2 bitmask code (`MAKE_OK`, in
https://github.com/ibethune/ap26/blob/master/src/cpu/cpuavx2.cpp --
PrimeGrid's AP27 search incrementally extends this same AP26 codebase).
That client ISA-specializes by duplicating the whole function per
architecture (a separate `cpusse41.cpp` reuses the same macro, but as a
second hand-adapted copy, not a shared core). The reference code anchors
its bitmasks to the *current* search position `A` and recomputes them live
as the search advances. This needs division at runtime, and projects
a bounded distance ahead per computation -- a reasonable design point
given its own constraints.

Our design indexes purely on residue: divisions are used only offline,
at table-build time, rather than repeatedly as the search moves. The
runtime hot path is unconditional lookup + AND, and on the CPU reference
that is the actual performance ceiling: bitmask-AND throughput. (The GPU
port's own bottleneck has not been separately profiled -- see "Merging
small primes" below.)

### Where each small prime is handled

Every prime up to 151 is covered, by one of four mechanisms:
- **2, 3:** gate `A mod 6` to 1 or 5. Since `6 | C`, every term is
  congruent to `A`.
- **29, 31:** also gated. When `p` does not divide `B`, the 27 terms
  occupy 27 distinct residues mod `p`, leaving `p-27` legal starts: 2
  mod 29, 4 mod 31, so `2*2*4 = 16` legal residues per `5394 = 6*29*31`.
  A batch (`512 = 32*16`) is exactly 32 whole periods, so every batch
  starts at the same residue mod 5394 and uses the same 512 offsets.
- **5..23** (the rest of `C`): every term is congruent to `A` mod `p`,
  so the only test is whether `p` divides `A`; merged into three
  tables mod 385, 299, 323.
- **37..151** (25 primes; "table primes" in the code): one table
  each. For candidate residue `r`, the term index hitting `p` is
  `i = -r * (B*C)^-1 mod p`; reject iff `i <= 26`.

This needs `B` coprime to 29, 31 and 37..151; otherwise a gate or table
degenerates, and both drivers refuse such a `B` at startup (about a
third of all `B` values).

### Why 32 primes

A larger table catches more composites per candidate, but every
prime checked costs a fixed lookup regardless of outcome. There is
no per-prime early-exit branching: the CPU path runs every lookup
unconditionally. The GPU path adds only one whole-invocation checkpoint
after its highest-yield primes, not a branch per prime -- see the GPU
port for why even that one checkpoint is worth it there. A prime `P`'s
marginal filtering value is approximately `27/P` for `P>27`---shrinking
fast as `P` grows---while its table-size cost grows linearly with `P`.
The practical optimum sits far below the "check with as many small
primes as you can" intuitive first choice: 32 primes (the 7 odd primes
`5..23` dividing `C`, plus the 25 table primes `37..151` -- `2, 3, 29,
31` are handled separately, by the gates above) already discard all but
7 of two million gated positions. On the CPU reference, the bottleneck
is AND throughput on a fixed-width bitmask, so limiting the number of
primes is beneficial there.

### Merging small primes: fewer lookups, slightly more memory used

The 7 odd primes `5..23` dividing `C` do not all need separate lookups.
A table indexed by the *product* of several coprime primes can encode
"bad mod any of them" in one entry (CRT), collapsing multiple lookups
into one at the cost of a bigger table (product-many bitmasks instead
of sum-many), while the tables' total size still remains modest. The
grouping used here---`{5,7,11}`, `{13,23}` and `{17,19}`---removes 4
lookups from a straightforward 7-bitmask implementation, while adding
+58KB of table entries.

With 28 lookups/512-bit batch---vs 32 without any
bitmask-merging---*~24% higher throughput* is sustained relative to
the 32-unmerged baseline on the CPU reference: more than the 32-vs-28
lookup-count ratio alone predicts (~14%, if every prime's lookup+update
cost the same). One plausible source of the extra ~10 points: the
grouped primes' index update is a real constant-modulus reduction
(`base_g[k] %= 385/299/323` in the code); if an unmerged baseline also
used a real modulus reduction per prime, merging would then remove four
of those costlier per-prime reductions, not just four lookups.

The same merge measures *~6.8% higher throughput* (444.7B -> 474.9B
candidates/sec) on the AMD Radeon 860M GPU port -- smaller than the
CPU's win. The GPU's own early-exit checkpoint (see the GPU port's
shader) covers the same primes in both forms, and its real measured
fire rate is ~82% (all 512 bits dead before the remaining lookups run);
a lookup-bound model from that fire rate predicts a ~27% win from
this merge, well above the measured 6.8%. So lookup+AND is not the
GPU's whole story -- something else (the remaining per-slot `% P`, the
per-invocation `atomicAdd(completed_count)`, or fixed per-invocation
overhead) likely dominates GPU time there; this has not been separately
profiled.

The same prime-product merge technique does *not* generalize to the 25
primes that don't divide `C` (37..151): table growth there scales with
the *product* of the grouped primes.

### Automatic SIMD-autovectorization

The choice of 512-bit vectors is not accidental: since the inner loop is
essentially ANDing 512-bit vectors, it does get recognized as amenable
to SIMD-ified processing on real ```gcc```/```clang``` builds -- but
some assembly required: the exact instructions depend on compiler,
version, and flags, so what follows is a reproducible sample of each,
not just a claim. `and_words_inplace`'s whole 16-word (512-bit) body
happens to exactly match one AVX-512 vector, so on real AVX-512
hardware the entire per-prime lookup+AND step (looping over the 25
table primes) collapses to a single `vpandd`:
```
gcc (SUSE Linux) 16.2.1, -O2 -march=x86-64-v4 -no-pie -g
(AMD Ryzen AI 7 PRO 350; x86-64-v4 is a portable AVX-512 ISA target,
not tied to this specific CPU)

		for (uint32_t p = 0; p < NB; p++) {
			const uint32_t *m = &rom_b[p][base_b[p] * WORDS];
  401520:	41 8b 0c 97          	mov    (%r15,%rdx,4),%ecx
	dst[  0 ] &= src[  0 ];
  401524:	48 8b 3c d5 20 51 40 	mov    0x405120(,%rdx,8),%rdi
  40152b:	00
		for (uint32_t p = 0; p < NB; p++) {
  40152c:	48 83 c2 01          	add    $0x1,%rdx
			const uint32_t *m = &rom_b[p][base_b[p] * WORDS];
  401530:	c1 e1 04             	shl    $0x4,%ecx
	dst[  0 ] &= src[  0 ];
  401533:	62 f1 7d 48 db 04 8f 	vpandd (%rdi,%rcx,4),%zmm0,%zmm0
		for (uint32_t p = 0; p < NB; p++) {
  40153a:	48 83 fa 19          	cmp    $0x19,%rdx
  40153e:	75 e0                	jne    401520
```
The shipped Makefile's default `-O2` (no `-march=`) does **not**
reproduce this: on the author's machine it emits 128-bit `pand` instead;
add `-march=x86-64-v4` to see the form above with gcc (Intel AVX-512
targets, and clang even at `x86-64-v4`, default to 256-bit vectors;
add `-mprefer-vector-width=512` to force). Even holding instruction
choice fixed, compiler tuning can still differ: the identical loop
built with `-march=znver5` (author's exact CPU) keeps the same `vpandd
zmm0` but picks `inc %rdx` over `add $0x1,%rdx` for the loop counter
-- same semantics, a different Zen5-specific scheduling preference.
Actual output will vary by compiler, version, and flags; this is only an
illustrative example.

The same source, unmodified, autovectorizes the identical loop into
NEON on ARM64 (Raspberry Pi 500, Cortex-A76) -- no intrinsics, no
per-architecture code path, confirmed on both compilers available
there:
```
gcc (Debian 12.2.0-14) 12.2.0, -O2

		for (uint32_t p = 0; p < NB; p++)
 36c:	d503201f 	nop
			and_words_inplace(combined, &(rom_b[p][base_b[p] * WORDS]));
 370:	b8404460 	ldr	w0, [x3], #4
 374:	f84084a2 	ldr	x2, [x5], #8
 378:	d37a6c00 	ubfiz	x0, x0, #6, #28
 37c:	8b000041 	add	x1, x2, x0
	dst[  0 ] &= src[  0 ];
 380:	3ce06844 	ldr	q4, [x2, x0]
 384:	ad409827 	ldp	q7, q6, [x1, #16]
 388:	3dc00c25 	ldr	q5, [x1, #48]
 38c:	4e241c63 	and	v3.16b, v3.16b, v4.16b
 390:	4e271c42 	and	v2.16b, v2.16b, v7.16b
 394:	4e261c21 	and	v1.16b, v1.16b, v6.16b
 398:	4e251c00 	and	v0.16b, v0.16b, v5.16b
		for (uint32_t p = 0; p < NB; p++)
 39c:	eb0300df 	cmp	x6, x3
 3a0:	54fffe81 	b.ne	370
```
Debian clang 14.0.6 at the same `-O2` produces the equivalent form
(both checked directly on this hardware, not assumed from one
compiler). NEON vectors are 128 bits, a quarter of AVX-512's 512, so
the identical 16-word loop body takes four `and v.16b` here where
x86-64 took one `vpandd` -- same portable source, same total work,
instruction count simply tracking native vector width.

## Layout

```
cpu/        Standalone CPU reference (no GPU needed). Start here.
            ap27_crt_sieve.c

gpu/        Vulkan compute port of the same design.
            ap27_rom28.comp        the compute shader
            ap27_rom28_kat.c       fixed-nbatches KAT verification driver
            ap27_rom28_search.c    real, unattended production search driver

presieve/   Second-stage filter (157..2677) for use between the sieve
            and an actual Fermat/Miller-Rabin primality test (not
            included in this repo), plus ap27_log_decimator -- a
            streaming filter for long unattended runs (see below).
            ap27_prefilter256.c    the filter itself, plus its self-test
            prefilter256_cli.c     standalone CLI front end
            ap27_log_decimator.c   streaming log-volume reducer
            common-hitline.h       shared hit-line parsing/printing

common-modinv.h   shared modular-inverse helper, used by both cpu/ and gpu/
common-cliargs.h  shared KEY=value CLI argument parser (common tooling)
ap27sample.py     standalone Python reference for the bitmask construction
                  shown above -- dependency-free, meant to be read
                  alongside the real C, not instead of it
testvectors/      real (not synthetic) captured hit-stream data used by
                  `make check`'s decimator regression test
```

## Building and running

```
make            # CPU sieve + presieve, no GPU toolchain needed
make gpu        # also build the GPU port (needs glslc + libvulkan-dev)
make check      # run every self-test, including 4 known-answer checks
```

Needs gcc or clang (the `__int128` the sieve's table construction
(`common-modinv.h`) and the presieve filter rely on is a compiler
extension, not standard C, and isn't available on e.g. MSVC).

Each component also runs standalone:

```
./cpu/ap27_crt_sieve KAT=<1..4>                   # verify a known AP27
./gpu/ap27_rom28_kat WINDOW=<n> NBATCHES=<n> KAT=<1..4>
echo "B=155368778 A=605185576317848261" | ./presieve/prefilter256_cli
```

### Known-answer tests

Every component is checked against four real, independently verifiable
AP27 examples -- not synthetic data:

| # | A | B |
|---|---|---|
| 1 | 605185576317848261 | 155368778 |
| 2 | 277699295941594831 | 170826477 |
| 3 | 224584605939537911 | 81292139 |
| 4 | 973044677408627339 | 174505864 |

A correct build finds each one at the exact predicted batch and bit
position, not just "eventually somewhere" -- that precision is itself
part of the check.

## Running a search

```
./cpu/ap27_crt_sieve --timed B=<b> A=<a_lo> SECONDS=<secs>
./gpu/ap27_rom28_search B=<b> A=<a_lo> SECONDS=<secs>
```

Every argument is a `KEY=value` token (shared parser:
`common-cliargs.h`), in any order. Each tool echoes its parsed values to
stderr, so a logged invocation line is self-documenting with nothing to
decode later. The same convention, via the same shared header, is used
by every multi-argument tool in this repository.

Both round `A_LO` up to the next batch start and scan forward. Each
32-prime survivor is one line, with `A = A0 + batch*172608 + offset[pos]`:

```
GPU: HIT B=<b> A=<decimal> status=unchecked batch=<n> pos=<s>
CPU: HIT B=..<last 5 digits> A=x<hex> ... id=<stream>
     (the full B is on the stream's first line)
```

`status=unchecked` means it survived the sieve only; it is not a
primality claim. The GPU driver reports progress, including the
current `A` for resuming, on stderr. Pipe either stream through
`./presieve/ap27_log_decimator [SAMPLE_EVERY=<n>]` before compressing --
it applies the second-stage filter inline, so only genuinely rare
survivors get logged. The second-stage filter's own survival ratio is
~1 in 193,060; the decimator also passes through a small sample of
ordinary hits, any new widest-arrival-gap record, and periodic
PROGRESS lines, so the real, observed reduction in archived bytes is
smaller than that and differs by driver: roughly 66,000x for GPU-driver
streams at default sampling, and roughly 1,800x for CPU-driver streams
(which emit PROGRESS far more often) -- both still large enough that
nothing a real analysis would need gets lost. "Widest-arrival-gap" is measured against a running maximum
`A`, not simple consecutive order, since the GPU's own parallel
invocations can report hits slightly out of order within a chunk; this
makes the record a watch on the arrival process itself (did hits stop
coming in the expected volume/shape), not a precise measurement of true
consecutive-`A` spacing:

```
<cpu or gpu search driver> | ./presieve/ap27_log_decimator [SAMPLE_EVERY=<n>] | zstd ...
```

To audit a build end to end, pipe live sieve output through
`./presieve/ap27_log_decimator --verify`: it re-checks every hit
against all primes up to 151 and exits nonzero on the first one that
should have been rejected (off by default; it costs extra CPU).

On a machine with several GPUs, set `AP27_GPU_INDEX=<n>` to pick the
Vulkan device (default: the first discrete GPU, else the first
integrated one).

## Results

Measured on real hardware (AMD Radeon 860M), the sieve alone -- no
presieve, no Fermat test -- sustains on the order of 10<sup>11</sup>
"candidates" (gated positions -- 512 per batch, not every integer
`A`) per second. The same source, unmodified, also runs correctly
on a second, independent GPU architecture and on an ARM64/NEON CPU
(Raspberry Pi 500's V3D GPU and Cortex-A76 CPU respectively) --
this confirms portability of the whole design. The 860M is an 8-CU
*integrated* GPU, not a discrete card -- if anything that makes the
number more notable, not less, since integrated GPUs are typically more
memory-bandwidth-constrained than discrete ones; the GPU port's own
actual bottleneck has not been separately profiled, though (see "Merging
small primes" above).

A production run scanned 1.45x10<sup>14</sup> candidates and its
measured survival rate matched the theoretical prediction to four
significant figures. At that candidate count, Poisson noise alone is
far too small to explain a mismatch of even one part in a thousand, so
a match this close is a real correctness signal: any genuine bug in the
rejection logic (a wrong table, a dropped prime, a sign error) would
show up clearly against noise this tight, not get lost in it.

None of this proves anything about the *density* of real AP27s in
unsearched territory -- that's a separate, much harder question this
repository doesn't address.

## License

BSD 3-Clause. See `LICENSE`.
