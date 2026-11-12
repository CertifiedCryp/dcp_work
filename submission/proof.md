# Unit-frame interval sieve for DCP over every Z_N: T = O(n^4 2^{(9/2)√n})

**Algorithm/configuration name:** `UFIS-9/2` (unit-frame interval sieve, b = ⌈2√n⌉, constants below).
All logarithms are base 2. N ≥ 2, n = ⌈log N⌉, s ∈ Z_N arbitrary, ω = e^{2πi/N}.
For a real x, frac(x) ∈ [0,1) and ‖x‖ = min_{z∈Z}|x − z| (circular distance, "turns").

**Credit.** The 2^{O(√n)} collision-sieve idea is Kuperberg's (SIAM J. Comput. 35(1), 2005); phase
reconstruction is Kitaev's (quant-ph/9511026). What this document adds is a self-contained
algorithm for *every* N, with exact uniformity in the sieve, no arbitrary-angle rotations
(hence no gate-synthesis error), worst-case resource bounds, and an explicit rational exponent
constant 9/2 under `dcp-uniform-clifford-t-bit-operations-v1`. We claim no improvement over the best known
DCP algorithms; the claim is a proof in this challenge's model.

## 1. Parameters (all computed classically in poly(n) time)

| symbol | definition |
|---|---|
| r | ⌊log N⌋, so 2^r ≤ N < 2^{r+1} |
| g | least integer ≥ 2 with gcd(g, N) = 1 (Lemma 1: g ≤ 4nL, L = ⌈log(n+1)⌉ + 1) |
| p | ⌈log(g+1)⌉, so g + 1 ≤ 2^p ≤ 2g + 2 |
| m | least integer with g^{m−1} ≥ N; m ≤ n + 1 |
| b | ⌈2√n⌉ (least integer with b² ≥ 4n) |
| d | ⌈(r − p − 1)/b⌉ if r > p + 1, else 0;  d ≤ n/b + 1 ≤ √n/2 + 1 |
| t | 32⌈log(64 m (d + p + 2))⌉ |
| M_d | max(2^{p+3} t, 2^b);  M_ℓ = 32^{d−ℓ} M_d for 0 ≤ ℓ ≤ d |
| Q | 4 m M_0 if d ≥ 1;  m (p+1) 2^{p+5} t if d = 0 |
| t_F | 32⌈log(64 (n+1) Q)⌉ |
| K | ⌈log(16 Q)⌉ + 1 |

Every quantity is a fixed function of N; no step depends on s except through measurement outcomes.

**Lemma 1 (small coprime base).** g ≤ 4nL with L = ⌈log(n+1)⌉ + 1.
*Proof.* For m' ≥ 1, I = ∫₀¹ x^{m'}(1−x)^{m'} dx = Σ_j c_j/(m'+j+1) with integers c_j, so
I·lcm(1..2m'+1) is a positive integer; since x(1−x) ≤ 1/4, I ≤ 4^{−m'} and lcm(1..2m'+1) ≥ 4^{m'}.
Hence lcm(1..x) ≥ 2^{x−2} for every integer x ≥ 1. Suppose every prime q ≤ x divides N. Then the number
of primes ≤ x is at most the number of distinct prime factors of N, at most log N ≤ n; and
lcm(1..x) = Π_{q≤x} q^{e_q} with q^{e_q} ≤ x, so lcm(1..x) ≤ x^n. With x = 4nL: log x = 2 + log n + log L ≤
2 + (L−1) + (L−1) = 2L, so x^n ≤ 2^{2nL} < 2^{4nL−2} ≤ lcm(1..x), a contradiction. So some prime ≤ 4nL
does not divide N. ∎  Consequently p = O(log n), 2^p = O(n log n), and m − 1 = ⌈log_g N⌉ ≤ n.

## 2. Primitives

**2.1 Fourier states.** Let U|y⟩ = |y+1 mod N⟩ on Z_N and |l̂⟩ = N^{−1/2} Σ_y ω^{−ly}|y⟩, so U|l̂⟩ = ω^l|l̂⟩ and
|x⟩ = N^{−1/2} Σ_l ω^{lx}|l̂⟩. With |ψ_l⟩ = (|0⟩ + ω^{ls}|1⟩)/√2, each source term is
(|0,x⟩ + |1,x+s⟩)/√2 = N^{−1/2} Σ_l ω^{lx} |ψ_l⟩|l̂⟩. Averaging over uniform x removes all l ≠ l' cross terms:

  ρ_s = (1/N) Σ_{l∈Z_N} |ψ_l⟩⟨ψ_l| ⊗ |l̂⟩⟨l̂|.   (1)

So a copy is exactly a uniformly random **true label** l ∈ Z_N together with |ψ_l⟩|l̂⟩; copies are independent.

**2.2 Hadamard tests (exact Clifford+T).** For an eigenstate |e⟩ of a unitary V with eigenvalue e^{2πiφ}:
ancilla |0⟩, H, controlled-V, (for a "Y-test" apply S†), H, measure. Outcome 0 has probability
(1+cos 2πφ)/2 (X-test) or (1+sin 2πφ)/2 (Y-test), and |e⟩ is unchanged. We use V = U^{c} for classical
constants c = 2^j mod N: controlled addition of a constant modulo N on an n-bit register, an exact
reversible circuit of O(n) X/CNOT/Toffoli gates and O(n) ancillas (add c into n+1 bits, subtract N,
copy the sign to a flag, conditionally add N back, clear the flag by comparing with c), with
Toffoli = 7 T gates plus Cliffords. On a phase qubit |ψ_λ⟩ (the qubit, not a register) the X- and Y-tests are
the single-qubit measurements of X and Y (H or S†H then measure): outcome 0 has probability
(1+cos 2πθ)/2, (1+sin 2πθ)/2 with θ = frac(λs/N).

**2.3 Octant estimator.** From t X-tests and t Y-tests with outcome-0 counts c₊, s₊, put
C = 2c₊ − t, S = 2s₊ − t (integers). The octant of (C,S) is decided by sign(C), sign(S) and |C| ≤ |S|
(ties to a fixed side); output the octant midpoint α = (2o+1)/16. If |C/t − cos 2πθ| ≤ 1/4 and
|S/t − sin 2πθ| ≤ 1/4, the point (C,S)/t is within √2/4 of the unit vector at angle 2πθ, so its angle differs
from 2πθ by at most arcsin(√2/4) < 0.362 rad = 0.0576 turns, and ‖α − θ‖ ≤ 1/16 + 0.0576 < 0.121 =: ε₀.
By Hoeffding (±1 outcomes, deviation 1/4), each of the two conditions fails with probability ≤ 2e^{−t/32}.

**Lemma 2 (base-c reconstruction).** Let c ≥ 2, θ_j = frac(c^j φ) for 0 ≤ j < J, and α_j with ‖α_j − θ_j‖ ≤ ε < 1/(2(c+1)).
Set θ̃_{J−1} = α_{J−1} and, for j = J−2..0, θ̃_j = the candidate frac((θ̃_{j+1} + ρ)/c), ρ ∈ {0..c−1}, closest to α_j.
Then ‖θ̃_j − θ_j‖ ≤ ε c^{−(J−1−j)}.
*Proof.* Induct downward with e = ‖θ̃_{j+1} − θ_{j+1}‖ ≤ ε. θ_j is one of the frac((θ_{j+1}+ρ)/c); the matching
candidate is within e/c of θ_j, any other is ≥ (1−e)/c from θ_j. Distances to α_j: ≤ ε + e/c versus
≥ (1−e)/c − ε; the first is smaller because 2ε + 2e/c ≤ 2ε(1 + 1/c) < 1/c. ∎
All candidates are rationals with denominator 16·2^{p}·g^{m} (or 16·2^{n}), so the reconstruction is exact
integer arithmetic on O(n) -bit numbers.

**2.4 Fourier sampling of a copy.** For j = 0..n run t_F X-tests and t_F Y-tests with V = U^{2^j mod N} on the
copy's Z_N register, form octant estimates α_j of frac(2^j l/N), reconstruct with Lemma 2 (c = 2, J = n+1,
ε₀ < 1/6), and record l* = round(N θ̃₀) mod N. If all 2(n+1) Hoeffding events hold, ‖θ̃₀ − l/N‖ ≤ ε₀/2^n < 1/(2N),
so l* = l. Then discard the register (by (1) it is in |l̂⟩, unentangled from the qubit given l).
P(l* ≠ l) ≤ 4(n+1) e^{−t_F/32} ≤ 4(n+1)/(64(n+1)Q) = 1/(16Q) (use e^{−y} ≤ 2^{−y}).

**2.5 Combination.** For phase qubits |ψ_{k₁}⟩|ψ_{k₂}⟩: CNOT (control 1st), measure the 2nd. Outcome 0 leaves
|ψ_{k₁+k₂}⟩, outcome 1 leaves |ψ_{k₁−k₂}⟩ up to global phase; each outcome has probability exactly 1/2, independent
of everything else. X on |ψ_k⟩ gives |ψ_{−k}⟩ up to global phase.

**2.6 Exact Bernoulli coins.** A Bernoulli(a/q) coin (integers 0 ≤ a ≤ q, q ≥ 1, w = ⌈log q⌉ so q > 2^{w−1}) uses up to K
rounds; each round succeeds with probability q/2^w > 1/2 and
draws w fair bits (H on |0⟩, measure) forming U ∈ [0,2^w) and stops if U < q, returning [U < a]. If all K
rounds fail, the whole run **aborts** (outputs 0). Given that it does not abort, the coin is exactly Bernoulli(a/q);
P(abort) < 2^{−K} ≤ 1/(32Q) per coin. Here q = W − |v| ≤ W ≤ 2^r, so w ≤ n.

## 3. Algorithm UFIS-9/2

For each frame j = 0..m−1 (fresh copies per frame), with unit a_j = g^{−j} mod N:

*d = 0 (only for bounded n).* For each i = 0..p, acquire 2^{p+5} t copies, Fourier-sample, keep those with
a_j l* ≡ −2^i (mod N), and abort if fewer than 2t remain.

*d ≥ 1.* Acquire 4M₀ copies; Fourier-sample each; frame label u = a_j l* mod N ∈ [0,N). Keep copies with u < 2^r,
discard others immediately; abort if fewer than M₀ remain, otherwise keep the first M₀ (acquisition order).
The list now holds M₀ phase qubits with frame labels in I₀ = [0, 2^r). Set ρ = r.
Level ℓ = 1..d: let β = min(b, ρ − p − 1) and W = 2^{ρ−β}. Stably sort the records (frame label, qubit address)
by bucket ⌊(u − lo(I))/W⌋ (2^β buckets of width W). Within each bucket pair records 1–2, 3–4, ... in list order.
For a pair (u₁,u₂): combine (2.5); on outcome 0 discard; on outcome 1 the new record is v = u₁ − u₂ ∈ (−W, W)
(an integer, *not* reduced mod N). Discard unless v ∈ [−W/2, W/2); then keep it with a coin of probability
(W/2)/(W − |v|) ≤ 1. Abort if fewer than M_ℓ records survive; keep the first M_ℓ. Set ρ ← ρ − β, I = [−W/2, W/2).
After level d, ρ = p + 1 and I_d = [−2^p, 2^p). For each i = 0..p, take records with v = −2^i; abort if fewer than 2t.
Apply X to the selected qubits.

*Estimation.* For each (j,i) use t selected qubits for X-tests and t for Y-tests; octant estimate α_{j,i}.
Lemma 2 (c = 2, J = p+1) gives β_j ≈ frac(g^j s/N); Lemma 2 (c = g, J = m) on (β_j) gives θ̃₀; output round(N θ̃₀) mod N.

**Label bookkeeping.** True label k and frame label u satisfy a_j k ≡ u (mod N) at all times: initially by
definition (when l* = l), and combination maps (k₁,k₂) to k₁ − k₂ with u₁ − u₂ ≡ a_j(k₁−k₂). A final record v = −2^i
has k ≡ −2^i g^j, so after X the qubit is |ψ_λ⟩ with λ = 2^i g^j and phase θ = frac(2^i g^j s/N).

## 4. Correctness

Define the **ideal process**: identical, except recorded labels equal true labels and coins never abort.
Couple real and ideal on the same randomness; they coincide unless some copy has l* ≠ l (prob ≤ Q/(16Q) = 1/16 by 2.4
and a union bound over Q copies) or some coin aborts (≤ 1/32, at most Q coins). Coin count: each level-ℓ pair uses
≤ 1 coin, and pairs ≤ list sizes ≤ 4M₀ per frame, total ≤ Q.

**Lemma 3 (uniformity).** In the ideal process, conditioned on everything observed about bucket indices, pairing,
combination outcomes and the surviving count, the frame labels of the list entering level ℓ are i.i.d. uniform
on an interval I of length 2^ρ.
*Proof.* Level 0: true labels are i.i.d. uniform on Z_N by (1); u = a_j l is a bijective image, hence uniform on [0,N);
conditioned on u < 2^r it is uniform on [0,2^r), and the selection (prefix by acquisition order) is independent of values.
Step: given i.i.d. uniform labels on I, condition on the bucket of each; offsets within buckets are i.i.d. uniform on
[0,W). Pairing depends only on buckets and list order. For a pair, v = u₁ − u₂ has P(v) = (W−|v|)/W² on (−W,W).
For v ∈ [−W/2, W/2) the coin gives P(v, kept) = (W−|v|)/W² · (W/2)/(W−|v|) = 1/(2W), constant over the W values,
so a kept v is uniform on [−W/2,W/2) (and kept with probability exactly 1/2). Pairs are disjoint, so kept values are
independent; their positions in the output list depend only on bucket indices and pair order, not on their values;
taking a prefix preserves i.i.d. uniformity. ∎

**Lemma 4 (no abort, ideal process).** Each of the following fails with probability ≤ e^{−t/2}:
(a) initial filter: X ~ Bin(4M₀, 2^r/N), 2^r/N > 1/2, mean > 2M₀; P(X < M₀) ≤ e^{−M₀/4} (Chernoff, δ = 1/2).
(b) level ℓ (input size M_{ℓ−1} ≥ 32 M_d ≥ 2^{b+1} ≥ 2·2^β): the number of pairs is P ≥ (M_{ℓ−1} − 2^β)/2 ≥ M_{ℓ−1}/4;
each pair survives independently with probability (1/2)(1/2) = 1/4; P(survivors < P/8) ≤ e^{−P/32} ≤ e^{−M_{ℓ−1}/128},
and P/8 ≥ M_{ℓ−1}/32 = M_ℓ.
(c) final: each record equals −2^i with probability 2^{−(p+1)} (Lemma 3), mean ≥ M_d 2^{−(p+1)} ≥ 4t;
P(< 2t) ≤ e^{−t/2}. Similarly in the case d = 0: probability 1/N > 2^{−(p+2)} per copy, mean ≥ 8t.
Since M₀, M_{ℓ−1}/128 ≥ M_d/4 ≥ 2^{p+1} t ≥ t/2, all bounds are ≤ e^{−t/2}. Per frame at most d + p + 2 events, so
P(any abort) ≤ m(d+p+2) e^{−t/2} ≤ m(d+p+2)·(64 m(d+p+2))^{−16} < 1/64.

**Lemma 5 (estimation).** Given no abort in the ideal process, each selected qubit for (j,i) is exactly |ψ_λ⟩,
λ = 2^i g^j. Each octant estimate is within ε₀ of θ_{j,i} = frac(2^i g^j s/N) except with probability
≤ 4e^{−t/32} ≤ 4/(64 m(d+p+2)); there are m(p+1) estimates, so all hold except with probability < 1/16.
Then by Lemma 2 (c = 2, J = p+1, ε₀ < 1/6), ‖β_j − frac(g^j s/N)‖ ≤ ε₀ 2^{−p} ≤ ε₀/(g+1) < 1/(2(g+1)),
and by Lemma 2 (c = g, J = m), ‖θ̃₀ − s/N‖ ≤ (1/(2(g+1))) g^{−(m−1)} < 1/(2N), so the output is s.

**Theorem.** For every N ≥ 2 and s, UFIS-9/2 outputs s with probability ≥ 1 − 1/16 − 1/32 − 1/64 − 1/16 > 3/4 > 2/3.

## 5. Resources (worst case, every branch; no restarts)

All counts are fixed by N; aborts only stop early. Write Λ = ⌈log(M + 2)⌉ for the address cost, where M is the
peak number of allocated bits and qubits. For n large enough that d ≥ 1 (all but finitely many n, since d = 0
forces N < 2^{p+2} = O(n log n)), using b ≥ 2√n, d ≤ √n/2 + 1, 2^{p+3} t = O(n log² n) = o(2^{2√n}):

- M_d ≤ 2^{2√n+1} + 2^{p+3}t = O(2^{2√n}),  32^d ≤ 32·2^{5n/b} ≤ 32·2^{(5/2)√n},  so M₀ = O(2^{(9/2)√n}).
- **Q** = 4 m M₀ ≤ 4(n+1) M₀ = **O(n 2^{(9/2)√n})**.
- t_F = O(log(nQ)) = O(√n); Λ = O(log(n M₀)) = O(√n); t, K = O(√n).

**S_Q.** Per frame at most 4M₀ stored phase qubits (rejected copies discarded immediately; lists only shrink),
plus one copy's n-qubit register and O(n) ancillas during sampling: **S_Q = O(2^{(9/2)√n})**.
**S_C.** Records of O(n) bits (frame label, qubit address, sort key) for ≤ 4M₀ entries, sorting buffer of the same size,
O(m(p+1)) estimates of O(n) bits, constants 2^j mod N, counters: **S_C = O(n 2^{(9/2)√n})**.

**T.** (i) Source: Q(n+1). (ii) Fourier sampling per copy: 2(n+1)t_F Hadamard tests, each O(n) gates, each gate
charged O(Λ) for addressing: O(n · √n · n · √n) = O(n³); classical counting, octants, Lemma 2 and the final
multiplication by N: O(n²) bit operations, each memory access charged O(Λ): O(n^{5/2}). Total O(Q n³).
(iii) Sieve: per frame and level, a merge sort of ≤ 4M₀ records costs O(M₀ log M₀) record operations, each O(n) bits
with O(Λ) addressing, i.e. O(M₀ · √n · n · √n) = O(n² M₀); d = O(√n) levels and m frames: O(m n^{5/2} M₀) = O(n^{5/2} Q).
Combinations, coins (≤ K rounds of ≤ n+1 fair bits each, comparisons): O(Q · √n · n · √n) = O(n² Q).
(iv) Estimation: O(m(p+1)t) single-qubit measurements and O(n²) arithmetic: negligible.
(v) Parameter computation (gcd trials up to 4nL, powers, logs): poly(n).
Hence **T = O(n³ Q) = O(n⁴ 2^{(9/2)√n})**.

| resource | bound | structured envelope | index |
|---|---|---|---|
| T | O(n⁴ 2^{(9/2)√n}) | 2^{(9/2) n^{1/2}} n⁴ | 1/2 |
| S_C | O(n 2^{(9/2)√n}) | 2^{(9/2) n^{1/2}} n | 1/2 |
| S_Q | O(2^{(9/2)√n}) | 2^{(9/2) n^{1/2}} | 1/2 |
| Q | O(n 2^{(9/2)√n}) | 2^{(9/2) n^{1/2}} n | 1/2 |

Index: for B(n) = 2^{c√n} n^k with c = 9/2, log log B = log(c√n + k log n) = (1/2) log n + O(1), so g(B) = 1/2.

## 6. Error tolerance

ε(n) = 0: the guarantee assumes every acquired copy is an independent ideal coset state ρ_s. No malformed-sample
robustness is claimed for this configuration.

## 7. Sanity check (not part of the proof)

`sim/dcp_m2.py` samples every measurement from its exact distribution under these parameters and recovers s in all
trials for small N, including the d = 0 and d = 1 regimes; larger N are checked with reduced t, M_d.
