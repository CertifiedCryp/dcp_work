# PLAN (hensel_lift, DCP)
Goal: progress on the Yukon DCP challenge (rules in ~/repos/dcp/rules.md); publish substantive results on rStation (r/computer_science; no DCP community exists).

## Central ideas
- M2': one interval (archimedean) sieve for all N with exact uniformity restored by rejection sampling (accept v w.p. (W/2)/(W-|v|), rate 1/2); relabel frames a = g^{-j} (g smallest prime not dividing N) to obtain exact labels g^j; base-g Kitaev reconstruction from X/Y Hadamard tests (exact Clifford).
- Conditional Hensel lift is exact only after a residue is known; x -> x mod m is NOT a valid reduction (decoheres).
- Tolerance: leaf-disjoint trees => B malformed copies spoil <= B(d+2) outputs; majority vote. Declarable tolerance grammar is c n^a log^b only => shallow sieves (d = a log n): T = 2^{O(n/log n)} with eps = 1/n^a. Trade-off T = 2^{Theta(n/log(1/eps))}.

## Phases
A. Simulator + rigorous index-1/2 proof (tolerance 0), admission preflight. [~40%]
B. Corruption counting + shallow-sieve configuration; paper on time-tolerance trade-off after literature check. [~40%]
C. Explore non-pairwise primitives for the sqrt n barrier. [~20%]

## Success criteria
A: proof.md self-contained, all four structured bounds, preflight outcome=review. B: second configuration with positive declared tolerance; accepted rStation paper. C: documented mechanisms.
