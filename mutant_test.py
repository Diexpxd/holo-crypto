"""Self-attack: does Theta (K-point evaluation over F_q) separate the Conway / Kinoshita-Terasaka mutants?

Needs a build of holocrypto_bench with kMaxN >= 13 (the default build caps at 12 crossings and rejects the KT word):

    sed 's/constexpr int kMaxN = 12; /constexpr int kMaxN = 16; /' holocrypto_bench.cpp > holocrypto_bench_n16.cpp
    zig c++ -O3 -march=native -std=c++17 holocrypto_bench_n16.cpp -o holocrypto_bench_n16.exe

Usage:  python mutant_test.py [path/to/holocrypto_bench_n16.exe]
"""
import random
import re
import subprocess
import sys
from pathlib import Path

EXE = sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).with_name("holocrypto_bench_n16.exe"))
Q = 2**31 - 1
CONWAY = [2, 2, 2, 1, -3, -2, -2, 1, -2, 1, -3]                 # 11n34 (CRC Concise Encyclopedia of Mathematics)
KT = [1, 1, 1, 3, 3, 2, -3, -1, -1, 2, -1, -3, -2]              # 11n42 (same source)
MAXLEN = 16


def run(word, k=8, seed=None):
    cmd = [EXE, "--verify", " ".join(map(str, word)), "-k", str(k)] + (["--seed", str(seed)] if seed is not None else [])
    out = subprocess.run(cmd, capture_output=True, text=True).stdout
    th = [int(m) for m in re.findall(r"theta\[\d+\]=(\d+)", out)]
    h = re.search(r"hash=([0-9a-f]{64})", out)
    return th, (h.group(1) if h else out.strip() or "ERR")


def equivalent_variants(base, count, rng):
    """Random words with the same braid closure: conjugation, far commutation, braid relation, s s^-1 insertion/removal."""
    def step(w):
        n = len(w)
        for _ in range(50):
            m = rng.choice(["rot", "comm", "braid", "ins", "del"])
            if m == "rot" and n > 1:
                r = rng.randrange(1, n); return w[r:] + w[:r]
            if m == "comm" and n > 1:
                i = rng.randrange(n - 1)
                if abs(abs(w[i]) - abs(w[i + 1])) >= 2:
                    v = w[:]; v[i], v[i + 1] = v[i + 1], v[i]; return v
            if m == "braid" and n > 2:
                i = rng.randrange(n - 2); a, b, c = w[i:i + 3]
                if a == c and abs(abs(a) - abs(b)) == 1 and (a > 0) == (b > 0):
                    return w[:i] + [b, a, b] + w[i + 3:]
            if m == "ins" and n + 2 <= MAXLEN:
                g = rng.choice([1, 2, 3]) * rng.choice([1, -1]); i = rng.randrange(n + 1)
                return w[:i] + [g, -g] + w[i:]
            if m == "del" and n > 1:
                i = rng.randrange(n - 1)
                if w[i] == -w[i + 1]:
                    return w[:i] + w[i + 2:]
        return w[:]
    seen, cur = {tuple(base)}, base[:]
    for _ in range(5000):
        if len(seen) >= count:
            break
        cur = step(cur)
        if len(cur) <= MAXLEN:
            seen.add(tuple(cur))
    return [list(x) for x in seen]


def main():
    print("== 1. Conway (11n34) vs Kinoshita-Terasaka (11n42)")
    for label, k, seed in (("K=8, default seed", 8, None), ("K=8, seed 12345", 8, 12345),
                           ("K=8, seed 0xDEADBEEF", 8, 0xDEADBEEF), ("K=16, default seed", 16, None)):
        c, hc = run(CONWAY, k, seed)
        t, ht = run(KT, k, seed)
        cm, _ = run([-x for x in CONWAY], k, seed)
        eq = sum(a == b for a, b in zip(c, t))
        neg = sum((a + b) % Q == 0 for a, b in zip(c, t))
        mirror_ok = all((a + b) % Q == 0 for a, b in zip(c, cm))
        print(f"  {label}: points={len(c)} | equal theta={eq} | opposite theta={neg} | "
              f"mirror control theta(K*)=-theta(K): {mirror_ok} | digests differ: {hc != ht}")
    print("\n== 2. Representation invariance (same knot, different braid words)")
    rng = random.Random(20261007)
    for name, base in (("Conway", CONWAY), ("Kinoshita-Terasaka", KT)):
        ref = run(base)[1]
        vs = equivalent_variants(base, 60, rng)
        bad = sum(run(v)[1] != ref for v in vs)
        print(f"  {name}: {len(vs)} equivalent words, digest mismatches: {bad}")


if __name__ == "__main__":
    main()
