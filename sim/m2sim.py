"""Classical simulation of M2' (rejection-uniform interval sieve + base-g frames + base-g Kitaev reconstruction).
A phase qubit is (label k, ideal state |0>+e^{2pi i k s/N}|1>); every measurement is sampled with its exact
quantum probability: CNOT+measure on (k1,k2) yields k1+k2 or k1-k2 w.p. 1/2 each; X/Y Hadamard tests give
+1 w.p. (1+cos)/2, (1+sin)/2."""
import math, numpy as np

def smallest_prime_not_dividing(N):
    g = 2
    while True:
        if all(g % p for p in range(2, int(g**0.5) + 1)) and N % g: return g
        g += 1

def sieve_frame(N, a, b, L0, rng):
    """L0 fresh coset samples -> number of phase qubits with frame label exactly 1, i.e. original label a^{-1}."""
    r = N.bit_length() - 1                     # R0 = 2^r <= N
    k = rng.integers(0, N, size=L0, dtype=np.int64)
    kp = (a * k) % N
    lab = kp[kp < (1 << r)]                    # uniform on [0, 2^r)
    lo, R, rem = 0, 1 << r, r - 1
    while R > 2:
        bb = min(b, rem); rem -= bb
        W = R >> bb
        bucket = (lab - lo) // W
        order = np.lexsort((rng.random(lab.size), bucket))
        lab, bucket = lab[order], bucket[order]
        # pair consecutive elements in the same bucket
        start = np.r_[0, np.flatnonzero(np.diff(bucket)) + 1]
        pos = np.arange(lab.size) - np.repeat(start, np.diff(np.r_[start, lab.size]))
        first = np.flatnonzero((pos % 2 == 0))
        first = first[first + 1 < lab.size]
        first = first[bucket[first] == bucket[first + 1]]
        v = lab[first] - lab[first + 1]
        keep = rng.random(v.size) < 0.5                         # CNOT outcome "difference"
        keep &= (v >= -(W // 2)) & (v < W // 2)
        keep &= rng.random(v.size) < (W / 2) / (W - np.abs(v))  # rejection -> uniform on [-W/2, W/2)
        lab = v[keep]; lo = -(W // 2); R = W
    return int(np.count_nonzero(lab == -1))                 # X-flip turns -1 into label 1

def estimate_theta(N, s, label, copies, rng):
    th = 2 * math.pi * ((label * s) % N) / N
    c = np.where(rng.random(copies) < (1 + math.cos(th)) / 2, 1, -1).mean()
    sn = np.where(rng.random(copies) < (1 + math.sin(th)) / 2, 1, -1).mean()
    return (math.atan2(sn, c) / (2 * math.pi)) % 1.0

def circ(x, y):
    d = (x - y) % 1.0
    return min(d, 1 - d)

def recover(N, s, b, rng, copies, L0):
    g = smallest_prime_not_dividing(N)
    m = math.ceil(math.log(N, g)) + 2
    used = 0; alphas = []
    for j in range(m):
        a = pow(g, -j, N)
        got = 0
        while got < 2 * copies:
            got += sieve_frame(N, a, b, L0, rng); used += L0
        alphas.append(estimate_theta(N, s, pow(g, j, N), copies, rng))
    theta = alphas[-1]
    for j in range(m - 2, -1, -1):
        theta = min((((theta + t) / g) % 1.0 for t in range(g)), key=lambda c: circ(c, alphas[j]))
    return round(N * theta) % N, used, g

if __name__ == '__main__':
    rng = np.random.default_rng(1)
    for N, b in [(64, 3), (97, 3), (100, 3), (1001, 4), (4096, 4), (65521, 5), (2**20, 5), (999983, 5)]:
        ok = used = 0; T = 20
        for _ in range(T):
            s = int(rng.integers(N))
            est, u, g = recover(N, s, b, rng, copies=80, L0=200000)
            ok += est == s; used += u
        print(f"N={N} g={g} b={b} success {ok}/{T} avg samples {used // T}", flush=True)
