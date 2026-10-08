"""Independent Python reference for the multipoint v2 digest; must match holocrypto_bench.exe --verify."""
import hashlib, sys
from holocrypto_engine_v1 import braid_to_upright, theta_eval_np

Q = 2**31 - 1
SEED = 0x484F4C4F43525950

def eval_points(K=8, seed=SEED):
    pts = []
    for k in range(K):
        h = hashlib.sha256(b"holocrypto-v2/point" + seed.to_bytes(8, "big") + bytes([k])).digest()
        pts.append((2 + int.from_bytes(h[:8], "big") % (Q - 3), 2 + int.from_bytes(h[8:16], "big") % (Q - 3)))
    return pts

def digest(word, m=4, K=8, seed=SEED):
    d = braid_to_upright(word, m)
    th = [int(theta_eval_np(d, a, b, Q)) for a, b in eval_points(K, seed)]
    msg = (b"holocrypto-v2/theta8" + Q.to_bytes(8, "big") + K.to_bytes(8, "big") + seed.to_bytes(8, "big")
           + b"".join(t.to_bytes(8, "big") for t in th))
    return th, hashlib.sha256(msg).hexdigest()

if __name__ == "__main__":
    w = [int(x) for x in sys.argv[1].split()]
    th, h = digest(w)
    for i, t in enumerate(th): print(f"theta[{i}]={t}")
    print("hash=" + h)
