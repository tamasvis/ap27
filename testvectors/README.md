# ap27_log_decimator real-data test vectors

`hitfilter_kat_survivor1.txt` -- 61 real, contiguous first-stage (32-prime
sieve) hit lines from an actual GPU run (`gpu/ap27_rom28_search`,
B=234222983), a real window of 30 hits before + 30 after one confirmed
second-stage-sieve survivor (line 31, A=288235970664735451). Verified
two ways before being saved: `./presieve/ap27_log_decimator` on this
exact file reports `survived=1` with that A, and independently, in
Python, that A has no term (A + i*234222983*223092870) divisible by any
prime up to 151 -- a real positive control, not synthetic data.

Expected regression-test result, any time ap27_log_decimator.c or
ap27_prefilter256.c change:

    ./presieve/ap27_log_decimator < testvectors/hitfilter_kat_survivor1.txt

Must report exactly `survived=1` in the FINAL summary line, with that
one survivor's A equal to 288235970664735451 (plain decimal -- this
tool's own output format; this file's own lines are the GPU driver's
decimal shape, see ap27_log_decimator.c's header comment for why the
CPU driver's compact hex lines are resolved differently). Any other
survived count is a real regression, not noise -- this window is small
and fixed, not subject to run-to-run Poisson variance the way a live
run is.
