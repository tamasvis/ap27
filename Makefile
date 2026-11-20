# SPDX-License-Identifier: BSD-3-Clause
#
# AP27 CRT/ROM-batch sieve: cpu/ is the standalone reference (build this
# first), gpu/ is a Vulkan port of the exact same design for wide
# parallel dispatch, and presieve/ is an optional second-stage filter
# that sits between either driver's hit stream and a real primality
# test. See the top-level README for why the design is structured this
# way; this file only lists how to build and check each piece.
#
# make            -- build the CPU sieve reference and the presieve
#                     (no Vulkan/GPU toolchain required)
# make gpu        -- also build the GPU port: ap27_rom28_kat (fixed-
#                     nbatches KAT verification) and ap27_rom28_search
#                     (real, unattended, duration-bound production search)
#                     (needs glslc + libvulkan-dev)
# make check      -- run every self-test/KAT check (CPU sieve, presieve,
#                     and GPU port if built). A GPU step that can't find a
#                     usable Vulkan driver/device prints SKIP and does not
#                     fail the overall result, since "no GPU on this
#                     machine" is a legitimate, common case. On a machine
#                     that's actually meant to exercise the GPU path (e.g.
#                     the real hardware used to validate this port), pass
#                     AP27_REQUIRE_GPU=1 to turn that same SKIP into a
#                     hard failure instead, so a real regression there
#                     (not just a missing GPU) can't go unnoticed as a
#                     quietly-passing SKIP.
# make clean

CC      ?= cc
CFLAGS  ?= -O2 -Wall -Wextra
GLSLC   ?= glslc
AP27_REQUIRE_GPU ?= 0

.PHONY: all gpu check clean

all: cpu/ap27_crt_sieve presieve/ap27_prefilter256_selftest presieve/prefilter256_cli \
	presieve/ap27_log_decimator

cpu/ap27_crt_sieve: cpu/ap27_crt_sieve.c common-modinv.h common-cliargs.h
	$(CC) $(CFLAGS) -o $@ cpu/ap27_crt_sieve.c

presieve/ap27_prefilter256_selftest: presieve/ap27_prefilter256.c
	$(CC) $(CFLAGS) -o $@ $<

# PREFILTER_MARCH exists because `-march=native` is gcc/x86-64 spelling
# only: clang doesn't support it at all (any target), and on ARM64 the
# equivalent is `-mcpu=native`. We don't consider this nice -- a single
# env var covering one flag isn't a real cross-compiler/cross-arch
# story -- but it's the only switch in this whole build that needs one,
# so it's good enough for a demonstration project: override it for your
# own compiler/arch (e.g. `PREFILTER_MARCH=-mcpu=native` for clang on
# ARM64, or empty it out entirely) rather than respecifying all of
# PREFILTER_CFLAGS just to drop this one flag.
PREFILTER_MARCH  ?= -march=native
PREFILTER_CFLAGS ?= -O3 $(PREFILTER_MARCH)

presieve/prefilter256_cli: presieve/prefilter256_cli.c presieve/ap27_prefilter256.c presieve/common-hitline.h
	$(CC) $(CFLAGS) $(PREFILTER_CFLAGS) -o $@ presieve/prefilter256_cli.c

presieve/ap27_log_decimator: presieve/ap27_log_decimator.c presieve/ap27_prefilter256.c presieve/common-hitline.h common-cliargs.h
	$(CC) $(CFLAGS) -o $@ presieve/ap27_log_decimator.c

gpu: gpu/ap27_rom28.spv gpu/ap27_rom28_kat gpu/ap27_rom28_search

gpu/ap27_rom28.spv: gpu/ap27_rom28.comp
	$(GLSLC) -o $@ $<

gpu/ap27_rom28_kat: gpu/ap27_rom28_kat.c common-modinv.h common-cliargs.h gpu/ap27_rom28.spv
	$(CC) $(CFLAGS) -o $@ gpu/ap27_rom28_kat.c -lvulkan

gpu/ap27_rom28_search: gpu/ap27_rom28_search.c common-modinv.h common-cliargs.h gpu/ap27_rom28.spv
	$(CC) $(CFLAGS) -o $@ gpu/ap27_rom28_search.c -lvulkan

check: all
	@if command -v $(GLSLC) >/dev/null 2>&1; then $(MAKE) gpu; fi
	@if [ "$(AP27_REQUIRE_GPU)" != "0" ] && { [ ! -x gpu/ap27_rom28_kat ] || [ ! -x gpu/ap27_rom28_search ]; }; then \
		echo "FAIL: AP27_REQUIRE_GPU is set but the GPU port is not built (glslc missing from PATH?)"; exit 1; fi
	@echo "=== CPU sieve: all 4 KATs ==="
	@for k in 1 2 3 4; do ./cpu/ap27_crt_sieve KAT=$$k || exit 1; done
	@echo "=== KEY=value parser: argument order is actually irrelevant, not just claimed (all 6 permutations each) ==="
	@ref=$$(./cpu/ap27_crt_sieve --timed B=234222983 A=288231422337521981 SECONDS=0 2>&1 >/dev/null | grep '^ap27_crt_sieve --timed:'); \
		[ -n "$$ref" ] || { echo "FAIL: reference ap27_crt_sieve --timed invocation produced no echo"; exit 1; }; \
		ok=1; \
		for perm in \
			'B=234222983 SECONDS=0 A=288231422337521981' \
			'A=288231422337521981 B=234222983 SECONDS=0' \
			'A=288231422337521981 SECONDS=0 B=234222983' \
			'SECONDS=0 B=234222983 A=288231422337521981' \
			'SECONDS=0 A=288231422337521981 B=234222983' \
		; do \
			out=$$(./cpu/ap27_crt_sieve --timed $$perm 2>&1 >/dev/null | grep '^ap27_crt_sieve --timed:'); \
			[ "$$out" = "$$ref" ] || { echo "FAIL: ap27_crt_sieve --timed order [$$perm] gave \"$$out\", expected \"$$ref\""; ok=0; }; \
		done; \
		[ $$ok -eq 1 ] && echo "PASS: ap27_crt_sieve --timed: all 6 orderings of B=/A=/SECONDS= parsed identically" || exit 1
	@if [ -x gpu/ap27_rom28_search ]; then \
		ref=$$(./gpu/ap27_rom28_search B=234222983 A=288231422337521981 SECONDS=0 2>&1 >/dev/null | grep '^ap27_rom28_search:'); \
		[ -n "$$ref" ] || { echo "FAIL: reference ap27_rom28_search invocation produced no echo"; exit 1; }; \
		ok=1; \
		for perm in \
			'B=234222983 SECONDS=0 A=288231422337521981' \
			'A=288231422337521981 B=234222983 SECONDS=0' \
			'A=288231422337521981 SECONDS=0 B=234222983' \
			'SECONDS=0 B=234222983 A=288231422337521981' \
			'SECONDS=0 A=288231422337521981 B=234222983' \
		; do \
			out=$$(./gpu/ap27_rom28_search $$perm 2>&1 >/dev/null | grep '^ap27_rom28_search:'); \
			[ "$$out" = "$$ref" ] || { echo "FAIL: ap27_rom28_search order [$$perm] gave \"$$out\", expected \"$$ref\""; ok=0; }; \
		done; \
		[ $$ok -eq 1 ] && echo "PASS: ap27_rom28_search: all 6 orderings of B=/A=/SECONDS= parsed identically" || exit 1; \
	fi
	@if [ -x gpu/ap27_rom28_kat ]; then \
		ref=$$(cd gpu && ./ap27_rom28_kat WINDOW=3000000 NBATCHES=20000 KAT=3 2>&1 >/dev/null | grep '^ap27_rom28_kat:'); \
		[ -n "$$ref" ] || { echo "FAIL: reference ap27_rom28_kat invocation produced no echo"; exit 1; }; \
		ok=1; \
		for perm in \
			'WINDOW=3000000 KAT=3 NBATCHES=20000' \
			'NBATCHES=20000 WINDOW=3000000 KAT=3' \
			'NBATCHES=20000 KAT=3 WINDOW=3000000' \
			'KAT=3 WINDOW=3000000 NBATCHES=20000' \
			'KAT=3 NBATCHES=20000 WINDOW=3000000' \
		; do \
			out=$$(cd gpu && ./ap27_rom28_kat $$perm 2>&1 >/dev/null | grep '^ap27_rom28_kat:'); \
			[ "$$out" = "$$ref" ] || { echo "FAIL: ap27_rom28_kat order [$$perm] gave \"$$out\", expected \"$$ref\""; ok=0; }; \
		done; \
		[ $$ok -eq 1 ] && echo "PASS: ap27_rom28_kat: all 6 orderings of WINDOW=/NBATCHES=/KAT= parsed identically" || exit 1; \
	fi
	@echo "=== presieve self-test ==="
	./presieve/ap27_prefilter256_selftest
	@echo "=== log decimator: real second-stage-survivor KAT ==="
	@./presieve/ap27_log_decimator < testvectors/hitfilter_kat_survivor1.txt 2>&1 >/dev/null \
		| grep -q "FINAL: total=61 survived=1 " \
		&& ./presieve/ap27_log_decimator < testvectors/hitfilter_kat_survivor1.txt 2>/dev/null \
			| grep -q "^HIT B=234222983 A=288235970664735451 status=unchecked stage2=survivor$$" \
		&& echo "PASS: decimator found exactly 1 survivor in the fixed 61-line window, with the expected A" \
		|| (echo "FAIL: log decimator KAT mismatch (count or survivor A)"; exit 1)
	@echo "=== log decimator: --verify on the same real window (upstream-sieve audit) ==="
	@./presieve/ap27_log_decimator --verify < testvectors/hitfilter_kat_survivor1.txt >/dev/null 2>&1 \
		&& echo "PASS: every line in the fixed window genuinely survives primes <= 151" \
		|| (echo "FAIL: --verify rejected a line in the fixed window -- see stderr"; \
			./presieve/ap27_log_decimator --verify < testvectors/hitfilter_kat_survivor1.txt >/dev/null; exit 1)
	@echo "=== CPU sieve: live output audited by --verify (--verify only checks a fixed file by itself; this closes the gap by running it against a real sieve) ==="
	@out=$$(mktemp); err=$$(mktemp); \
		./cpu/ap27_crt_sieve --timed B=234222983 A=288231422337521981 SECONDS=2 >"$$out"; src_rc=$$?; \
		if [ $$src_rc -ne 0 ]; then echo "FAIL: ap27_crt_sieve --timed exited $$src_rc"; rm -f "$$out" "$$err"; exit 1; fi; \
		./presieve/ap27_log_decimator --verify <"$$out" >/dev/null 2>"$$err"; rc=$$?; \
		if [ $$rc -eq 0 ]; then echo "PASS: every live CPU hit survives primes <= 151"; \
		else echo "FAIL"; cat "$$err"; rm -f "$$out" "$$err"; exit 1; fi; \
		rm -f "$$out" "$$err"
	@if [ -x gpu/ap27_rom28_search ]; then \
		echo "=== GPU search driver: live output audited by --verify ==="; \
		out=$$(mktemp); err=$$(mktemp); \
		./gpu/ap27_rom28_search B=234222983 A=288231422337521981 SECONDS=5 >"$$out"; drc=$$?; \
		if [ $$drc -eq 77 ] && [ "$(AP27_REQUIRE_GPU)" != "0" ]; then echo "FAIL: no usable GPU for ap27_rom28_search, and AP27_REQUIRE_GPU is set"; rm -f "$$out" "$$err"; exit 1; \
		elif [ $$drc -eq 77 ]; then echo "SKIP: no usable GPU for ap27_rom28_search"; rm -f "$$out" "$$err"; \
		elif [ $$drc -ne 0 ]; then echo "FAIL: ap27_rom28_search exited $$drc"; rm -f "$$out" "$$err"; exit 1; \
		else \
			./presieve/ap27_log_decimator --verify <"$$out" >/dev/null 2>"$$err"; rc=$$?; \
			if [ $$rc -eq 0 ]; then echo "PASS: every live GPU hit survives primes <= 151"; \
			else echo "FAIL"; cat "$$err"; rm -f "$$out" "$$err"; exit 1; fi; \
			rm -f "$$out" "$$err"; \
		fi; \
	fi
	@if [ -x gpu/ap27_rom28_kat ]; then \
		echo "=== GPU port: all 4 KATs ==="; \
		( cd gpu && for k in 1 2 3 4; do \
			./ap27_rom28_kat WINDOW=3000000 NBATCHES=20000 KAT=$$k; rc=$$?; \
			if [ $$rc -eq 77 ]; then exit 77; fi; \
			if [ $$rc -ne 0 ]; then exit 1; fi; \
		done ); rc=$$?; \
		if [ $$rc -eq 77 ] && [ "$(AP27_REQUIRE_GPU)" != "0" ]; then echo "FAIL: no usable GPU for ap27_rom28_kat, and AP27_REQUIRE_GPU is set"; exit 1; \
		elif [ $$rc -eq 77 ]; then echo "SKIP: no usable GPU for ap27_rom28_kat"; \
		elif [ $$rc -ne 0 ]; then exit 1; fi; \
	fi
	@if [ -x gpu/ap27_rom28_kat ]; then \
		echo "=== GPU port: all 4 KATs past a real chunk boundary ==="; \
		echo "(CHUNK_BATCHES=2500000; this run is 2600000 batches, so it genuinely"; \
		echo "exercises 2 chunks and the chunk-boundary reconstruction/overlap logic,"; \
		echo "not just the single-chunk path the KATs above stay within)"; \
		( cd gpu && for k in 1 2 3 4; do \
			./ap27_rom28_kat WINDOW=431692608000 NBATCHES=2600000 KAT=$$k; rc=$$?; \
			if [ $$rc -eq 77 ]; then exit 77; fi; \
			if [ $$rc -ne 0 ]; then exit 1; fi; \
		done ); rc=$$?; \
		if [ $$rc -eq 77 ] && [ "$(AP27_REQUIRE_GPU)" != "0" ]; then echo "FAIL: no usable GPU for ap27_rom28_kat, and AP27_REQUIRE_GPU is set"; exit 1; \
		elif [ $$rc -eq 77 ]; then echo "SKIP: no usable GPU for ap27_rom28_kat"; \
		elif [ $$rc -ne 0 ]; then exit 1; fi; \
	fi
	@echo "ALL CHECKS PASSED"

clean:
	rm -f cpu/ap27_crt_sieve presieve/ap27_prefilter256_selftest \
		presieve/prefilter256_cli presieve/ap27_log_decimator \
		gpu/ap27_rom28.spv gpu/ap27_rom28_kat gpu/ap27_rom28_search
