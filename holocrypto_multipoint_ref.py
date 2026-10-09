"""Independent Python reference for the multipoint v2 digest; must match holocrypto_bench.exe --verify."""
import hashlib, sys
from holocrypto_engine_v1 import braid_to_upright, theta_eval_np

Q = 2**31 - 1
SEED = 0x484F4C4F43525950

def check_params(K, seed):
    """K = 0 made digest() a CONSTANT (same value for every knot); K > 255 or a bad seed crashed obscurely."""
    if isinstance(K, bool) or not isinstance(K, int) or not 1 <= K <= 255:
        raise ValueError(f"K must be an int in [1, 255], got {K!r}")
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed < 2**64:
        raise ValueError(f"seed must be an int in [0, 2**64), got {seed!r}")


def eval_points(K=8, seed=SEED):
    check_params(K, seed)
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


class DegenerateKey(ValueError):
    """theta = 0 at every point: unknot, amphichiral knots, K#K*... all of them share ONE digest."""


def digest_checked(word, m=4, K=8, seed=SEED):
    """Same value as digest(), but refuses the degenerate theta = 0 class instead of returning the one
    constant digest shared by the unknot, the figure-eight knot, every K # K*, etc. (README, section 3.6)."""
    th, h = digest(word, m, K, seed)
    if not any(th):
        raise DegenerateKey("theta = 0 at all points: unknot / amphichiral / K#K*; not a usable key")
    return th, h


def digest_v3(word, m=4, K=8, seed=SEED):
    """OPTIONAL v3 digest (different domain tag, NOT compatible with digest()). It also hashes
    Delta(T1), Delta(T2), Delta(T3) at every point, so K # K* (theta = 0) no longer collapses onto
    the unknot digest. It does NOT fix colliding knots that share Delta (Conway / Kinoshita-Terasaka
    have Delta = 1) nor the connected-sum families built from them. Cost: ~0 extra (Delta is already
    computed inside theta)."""
    d = braid_to_upright(word, m)
    parts = []
    for a, b in eval_points(K, seed):
        th, dl = theta_eval_np(d, a, b, Q, return_deltas=True)
        parts += [*dl, int(th)]
    msg = (b"holocrypto-v3/alexander+theta" + Q.to_bytes(8, "big") + K.to_bytes(8, "big") + seed.to_bytes(8, "big")
           + b"".join(x.to_bytes(8, "big") for x in parts))
    return parts, hashlib.sha256(msg).hexdigest()


if __name__ == "__main__":
    w = [int(x) for x in sys.argv[1].split()]
    th, h = digest(w)
    for i, t in enumerate(th): print(f"theta[{i}]={t}")
    print("hash=" + h)
