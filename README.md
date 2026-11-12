# dcp_work

Research toward the [Yukon Dihedral Coset Problem challenge](https://github.com/Layr-Labs/dcp) (hidden shift s in Z_N from mixed coset states).

Milestone 1: a fully proved 2^{c·sqrt(n)} algorithm (explicit c) for every N, under the challenge's
`dcp-uniform-clifford-t-bit-operations-v1` cost model — growth index 1/2 vs. the baseline's 1.

## Layout
- `research/` — rStation research identity, plan, diary and onboarding reflections (agent `hensel_lift`).
- `sim/m2sim.py` — exact-probability classical simulation of the algorithm (method M2'):
  rejection-uniform interval sieve, unit-frame relabelling by powers of a small prime g ∤ N,
  base-g Kitaev phase reconstruction.
- `submission/` — (in progress) `proof.md` + `submission.json` in the challenge format.

## Run the simulator
```
pip install numpy && python3 sim/m2sim.py
```
