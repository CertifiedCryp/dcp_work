"""Exact-probability simulation of the submitted algorithm (method M2'), mirroring submission/proof.md.

Every quantum step is replaced by sampling from its exact outcome distribution:
  * Fourier sampling of a coset copy: label l uniform in Z_N; the Hadamard tests on controlled
    U^{2^j mod N} give +1 w.p. (1+cos 2pi phi)/2 (X) and (1+sin 2pi phi)/2 (Y), phi = frac(2^j l/N).
  * CNOT+measure on phase qubits (k1,k2): outcome 1 ("difference", label k1-k2) w.p. exactly 1/2.
  * X/Y tests on a phase qubit with label lam: phi = frac(lam s / N).
Labels are tracked both as recorded (estimated) and true values, so Fourier-sampling errors propagate honestly.
Parameters default to the proof's values; `scale` overrides shrink them for larger N (logic check only)."""
import math, numpy as np

def isqrt_ceil(n): r = math.isqrt(n); return r if r * r == n else r + 1
def clog2(x): return (x - 1).bit_length()             # ceil(log2 x), x >= 1

def params(N, t_override=None, md_override=None):
    n = clog2(N); r = N.bit_length() - 1
    g = 2
    while math.gcd(g, N) != 1: g += 1
    p = g.bit_length()                                 # ceil(log2(g+1))
    m = 1
    while g ** (m - 1) < N: m += 1
    b = 1
    while b * b < 4 * n: b += 1                       # b = ceil(2 sqrt n)
    d = -(-(r - p - 1) // b) if r > p + 1 else 0
    t = 32 * clog2(64 * m * (d + p + 2))
    if t_override: t = t_override
    Md = max(2 ** (p + 3) * t, 2 ** b)
    if md_override: Md = md_override
    M0 = 32 ** d * Md
    Q = m * 4 * M0 if d > 0 else m * (p + 1) * 2 ** (p + 5) * t
    tF = 32 * clog2(64 * (n + 1) * Q)
    return dict(n=n, r=r, g=g, p=p, m=m, b=b, d=d, t=t, Md=Md, M0=M0, Q=Q, tF=tF, K=clog2(16 * Q))

def octant(cp, sp, t):
    """cp, sp: counts of +1 outcomes among t X/Y tests -> midpoint of the octant of (c_hat, s_hat), in turns."""
    c = 2 * cp - t; s = 2 * sp - t                     # t * c_hat, t * s_hat (integers)
    ang_quadrant = np.where(s >= 0, np.where(c > 0, 0, 1), np.where(c <= 0, 2, 3))
    # within quadrant, first or second half by |c| vs |s| (rotated appropriately)
    ac, as_ = np.abs(c), np.abs(s)
    first = np.where((ang_quadrant % 2) == 0, as_ <= ac, ac <= as_)
    octn = 2 * ang_quadrant + np.where(first, 0, 1)
    return (2 * octn + 1) / 16.0

def circ(x, y):
    dd = np.mod(x - y, 1.0); return np.minimum(dd, 1 - dd)

def chain(alphas, c):
    """alphas[..., j] estimates frac(c^j phi); backward base-c reconstruction of phi."""
    th = alphas[..., -1]
    for j in range(alphas.shape[-1] - 2, -1, -1):
        cands = np.stack([np.mod((th + r) / c, 1.0) for r in range(c)], -1)
        idx = np.argmin(circ(cands, alphas[..., j:j + 1]), -1)
        th = np.take_along_axis(cands, idx[..., None], -1)[..., 0]
    return th

def tests(phi, t, rng):
    cp = rng.binomial(t, (1 + np.cos(2 * np.pi * phi)) / 2)
    sp = rng.binomial(t, (1 + np.sin(2 * np.pi * phi)) / 2)
    return cp, sp

def fourier_sample(N, cnt, P, rng, chunk=1 << 20):
    if cnt > chunk:
        parts = [fourier_sample(N, min(chunk, cnt - i), P, rng) for i in range(0, cnt, chunk)]
        return np.concatenate([a for a, _ in parts]), np.concatenate([b for _, b in parts])
    n = P['n']; l = rng.integers(0, N, size=cnt, dtype=np.int64)
    pw = np.array([pow(2, j, N) for j in range(n + 1)], dtype=np.int64)
    phi = np.mod((l[:, None] * pw[None, :]) % N / N, 1.0)
    alphas = octant(*tests(phi, P['tF'], rng), P['tF'])
    est = np.rint(N * chain(alphas, 2)).astype(np.int64) % N
    return est, l

def frame(N, s, j, P, rng, stats):
    """Returns dict i -> array of true labels of phase qubits whose recorded label is 2^i g^j (or None on abort)."""
    g, p, d, b, t = P['g'], P['p'], P['d'], P['b'], P['t']
    a = pow(g, -j, N)
    out = {}
    if d == 0:
        for i in range(p + 1):
            cnt = 2 ** (p + 5) * t; stats['Q'] += cnt
            rec, tru = fourier_sample(N, cnt, P, rng)
            sel = (a * rec) % N == (-(2 ** i)) % N
            if sel.sum() < 2 * t: return None
            out[i] = (-tru[sel])[:2 * t] % N                  # X-flip negates the label
        return out
    cnt = 4 * P['M0']; stats['Q'] += cnt
    rec, tru = fourier_sample(N, cnt, P, rng)
    stats['mislabel'] += int((rec != tru).sum())
    u = (a * rec) % N
    r = P['r']; keep = u < (1 << r)
    if keep.sum() < P['M0']: return None
    lab, tl = u[keep][:P['M0']], tru[keep][:P['M0']]
    lo, rho, M = 0, r, P['M0']
    while rho > p + 1:
        bb = min(b, rho - (p + 1)); W = 1 << (rho - bb)
        bucket = (lab - lo) // W
        order = np.argsort(bucket, kind='stable')               # within bucket: original list order
        lab, tl, bucket = lab[order], tl[order], bucket[order]
        start = np.r_[0, np.flatnonzero(np.diff(bucket)) + 1]
        pos = np.arange(lab.size) - np.repeat(start, np.diff(np.r_[start, lab.size]))
        f = np.flatnonzero(pos % 2 == 0); f = f[f + 1 < lab.size]; f = f[bucket[f] == bucket[f + 1]]
        v = lab[f] - lab[f + 1]; tv = (tl[f] - tl[f + 1]) % N
        ok = rng.random(v.size) < 0.5
        ok &= (v >= -(W // 2)) & (v < W // 2)
        ok &= rng.random(v.size) < (W / 2) / (W - np.abs(v))
        M //= 32
        if ok.sum() < M: return None
        lab, tl = v[ok][:M], tv[ok][:M]; lo, rho = -(W // 2), rho - bb
    for i in range(p + 1):
        sel = lab == -(2 ** i)
        if sel.sum() < 2 * t: return None
        out[i] = (-tl[sel])[:2 * t] % N
    return out

def run(N, s, P, rng):
    g, p, m, t = P['g'], P['p'], P['m'], P['t']
    stats = {'Q': 0, 'mislabel': 0}
    th_hat = []
    for j in range(m):
        out = frame(N, s, j, P, rng, stats)
        if out is None: return None, stats
        al = []
        for i in range(p + 1):
            tl = out[i]
            phi = np.mod((tl * s) % N / N, 1.0)
            cp = rng.binomial(1, (1 + np.cos(2 * np.pi * phi[:t])) / 2).sum()
            sp = rng.binomial(1, (1 + np.sin(2 * np.pi * phi[t:2 * t])) / 2).sum()
            al.append(octant(np.array(cp), np.array(sp), t))
        th_hat.append(chain(np.array(al, dtype=float), 2))
    th0 = chain(np.array(th_hat, dtype=float), g)
    return int(np.rint(N * th0)) % N, stats

if __name__ == '__main__':
    import sys
    rng = np.random.default_rng(7)
    print("-- exact proof parameters --")
    for N, T in [(2, 20), (3, 20), (5, 20), (6, 20), (12, 20), (30, 20), (64, 5), (210, 3)]:
        P = params(N); ok = 0; mis = 0
        for _ in range(T):
            s = int(rng.integers(N)); est, st = run(N, s, P, rng); ok += est == s; mis += st['mislabel']
        print(N, {k: P[k] for k in 'g p m b d t Md Q tF'.split()}, f"success {ok}/{T} mislabels {mis}", flush=True)
    print("-- reduced t, M_d (sieve/reconstruction logic at larger N; not the proof's constants) --")
    for N in [1001, 4096, 30030, 65521, 2**20, 999983, 510510, 2**24 - 3]:
        P = params(N, t_override=24, md_override=3 * 2 ** (P0 := params(N))['p'] * 24)
        P['tF'] = 64
        ok = 0; T = 10; fails = 0
        for _ in range(T):
            s = int(rng.integers(N))
            P2 = dict(P); P2['M0'] = int(P['Md'] * 9 ** P['d'])   # observed per-level yield ~1/7.5 > 1/9
            est, st = run(N, s, P2, rng); ok += est == s; fails += est is None
        print(N, {k: P[k] for k in 'g p m b d'.split()}, f"success {ok}/{T} aborts {fails}", flush=True)
