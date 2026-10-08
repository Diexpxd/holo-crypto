"""Fraction of braid words (closure = one knot) with theta = 0 at 2 random points of F_q, by (m strands, n crossings).
Words are uniform among those of length n with an m-cyclic permutation (sampling via _cyclic_word)."""
import random, sys, time, json
from holocrypto_engine_v1 import _cyclic_word, braid_to_upright, theta_nonzero

rng = random.Random(20261007)
grids = {4: [5, 7, 9, 11, 15, 21, 31, 51], 6: [5, 7, 9, 11, 15, 21, 31, 51], 10: [9, 11, 15, 21, 31, 51, 100 + 1]}
NSAMPLES = int(sys.argv[1]) if len(sys.argv) > 1 else 300
res = []
for m, ns in grids.items():
    for n in ns:
        t0 = time.time(); zeros = 0
        N = NSAMPLES if n <= 51 else max(100, NSAMPLES // 2)
        for _ in range(N):
            w = _cyclic_word(n, m, rng)
            if not theta_nonzero(braid_to_upright(w, m), rng):
                zeros += 1
        res.append(dict(m=m, n=n, N=N, zeros=zeros, frac=zeros / N, secs=round(time.time() - t0, 1)))
        print(f"m={m:2d} n={n:3d} N={N:3d} theta=0: {zeros:3d} ({100*zeros/N:5.1f} %)  [{time.time()-t0:.0f}s]", flush=True)
json.dump(res, open("theta_zero_density.json", "w"), indent=1)
