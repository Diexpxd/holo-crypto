// holocrypto_bench.cpp - C++17 port of holocrypto_engine_v1.py (theta_eval_np + holocrypto_hash)
// and a performance / combinatorics benchmark over 4-strand braids.
//
// Mathematics (identical to the .py): theta(D)(T1,T2) in F_q, q = 2^31-1,
//   A(T) = I - sum_c ( T^s E_{i,i+1} + (1-T^s) E_{i,j+1} + E_{j,j+1} ),  G = A^-1,  T3 = T1*T2
//   theta = Delta(T1) Delta(T2) Delta(T3) * (sum F1 + sum F2 + sum F3)
//
// v2 (multipoint): theta is evaluated at K points (T1_k, T2_k), k=0..K-1, derived from a public seed:
//   (T1_k, T2_k) = two 8-byte integers of SHA-256("holocrypto-v2/point" || seed || k), reduced to [2, q-2].
//   digest = SHA-256("holocrypto-v2/theta8" || q || K || seed || theta_0 || ... || theta_{K-1}),
//   8-byte big-endian integers. K=8 by default (8 x 31 = 248 bits of theta output space).
//   The single-point version (v1) is kept in holocrypto_bench_v1_single_point.cpp.
//
// Memory: everything lives on the stack with fixed sizes (no std::vector on the hot path):
// uint32_t matrices of kMaxM x kMaxM, products in uint64_t, Mersenne reduction without division.
//
// Usage: holocrypto_bench [-t THREADS] [-k POINTS] [--seed N] [--lens 5,6,7,8]
//        holocrypto_bench --verify "1 -2 3 1 2"        (prints the K thetas and the hash of a word)
//        holocrypto_bench --verify "1 -2 3 1 2" --v3   (optional v3 digest: also hashes Delta(T1..T3) per point)
//        holocrypto_bench --selftest                    (SHA-256 test vectors)
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace {

// ============================================================================
// 1. Arithmetic in F_q, q = 2^31 - 1 (Mersenne prime)
// ============================================================================
constexpr uint64_t Q = 0x7FFFFFFFULL;
using u32 = uint32_t;
using u64 = uint64_t;

inline u32 red(u64 x) {                       // x < 2^62  ->  x mod Q
    x = (x & Q) + (x >> 31);
    x = (x & Q) + (x >> 31);
    return static_cast<u32>(x >= Q ? x - Q : x);
}
inline u32 add(u32 a, u32 b) { u32 r = a + b; return r >= Q ? r - static_cast<u32>(Q) : r; }
inline u32 sub(u32 a, u32 b) { return a >= b ? a - b : a + static_cast<u32>(Q) - b; }
inline u32 neg(u32 a) { return a ? static_cast<u32>(Q) - a : 0; }
inline u32 mul(u32 a, u32 b) { return red(static_cast<u64>(a) * b); }
inline u32 fromInt(long long v) {
    long long r = v % static_cast<long long>(Q);
    return static_cast<u32>(r < 0 ? r + static_cast<long long>(Q) : r);
}
u32 power(u32 b, u64 e) {
    u32 r = 1;
    while (e) { if (e & 1) r = mul(r, b); b = mul(b, b); e >>= 1; }
    return r;
}
inline u32 inv(u32 a) { return power(a, Q - 2); }     // a != 0

// ============================================================================
// 2. SHA-256 (FIPS 180-4)
// ============================================================================
struct Sha256 {
    static constexpr u32 K[64] = {
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};

    static inline u32 rotr(u32 x, int n) { return (x >> n) | (x << (32 - n)); }

    static void block(u32 h[8], const uint8_t* p) {
        u32 w[64];
        for (int i = 0; i < 16; ++i)
            w[i] = (u32(p[4 * i]) << 24) | (u32(p[4 * i + 1]) << 16) | (u32(p[4 * i + 2]) << 8) | p[4 * i + 3];
        for (int i = 16; i < 64; ++i) {
            u32 s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
            u32 s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }
        u32 a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], hh = h[7];
        for (int i = 0; i < 64; ++i) {
            u32 t1 = hh + (rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)) + ((e & f) ^ (~e & g)) + K[i] + w[i];
            u32 t2 = (rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)) + ((a & b) ^ (a & c) ^ (b & c));
            hh = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
        }
        h[0] += a; h[1] += b; h[2] += c; h[3] += d; h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
    }

    static void hash(const uint8_t* msg, size_t len, uint8_t out[32]) {
        u32 h[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
        size_t full = len / 64;
        for (size_t i = 0; i < full; ++i) block(h, msg + 64 * i);
        uint8_t tail[128] = {0};
        size_t rem = len - 64 * full;
        std::memcpy(tail, msg + 64 * full, rem);
        tail[rem] = 0x80;
        size_t tlen = rem + 1 + 8 <= 64 ? 64 : 128;
        u64 bits = static_cast<u64>(len) * 8;
        for (int i = 0; i < 8; ++i) tail[tlen - 1 - i] = static_cast<uint8_t>(bits >> (8 * i));
        for (size_t o = 0; o < tlen; o += 64) block(h, tail + o);
        for (int i = 0; i < 8; ++i) {
            out[4 * i] = static_cast<uint8_t>(h[i] >> 24); out[4 * i + 1] = static_cast<uint8_t>(h[i] >> 16);
            out[4 * i + 2] = static_cast<uint8_t>(h[i] >> 8); out[4 * i + 3] = static_cast<uint8_t>(h[i]);
        }
    }
};
constexpr u32 Sha256::K[64];

std::string hex(const uint8_t* d, size_t n) {
    static const char* H = "0123456789abcdef";
    std::string s;
    for (size_t i = 0; i < n; ++i) { s += H[d[i] >> 4]; s += H[d[i] & 15]; }
    return s;
}

// ============================================================================
// 3. Upright diagram (long knot) from a closed braid
// ============================================================================
constexpr int kMaxN = 12;                 // max crossings
constexpr int kMaxM = 2 * kMaxN + 1;      // max edges

struct Crossing { int s, i, j; };         // sign, over edge, under edge (from 1)

struct Knot {
    int n = 0;
    std::array<Crossing, kMaxN> x{};
    std::array<int, kMaxM> rot{};         // phi_1..phi_{2n+1}
    int expDelta() const {                // (-sum(phi) - w) / 2, integer division towards -inf like Python's divmod
        int sphi = 0, w = 0;
        for (int k = 0; k < 2 * n + 1; ++k) sphi += rot[k];
        for (int c = 0; c < n; ++c) w += x[c].s;
        int v = -sphi - w;
        return (v >= 0) ? v / 2 : -((-v + 1) / 2);
    }
};

// Closure of the braid (letters +-k = sigma_k^{+-1}, m strands) cut at the closing arc of strand 1.
// Mirrors _walk_braid + braid_to_upright. Returns false if the closure is not a knot (several components).
bool braidToKnot(const int8_t* word, int n, int m, Knot& K) {
    if (n > kMaxN) return false;
    int stepC[2 * kMaxN], stepOver[2 * kMaxN], stepArcs[2 * kMaxN];
    int steps = 0, arcs = 0, strands = 0, x = 1;
    for (;;) {
        int pos = x;
        for (int c = 0; c < n; ++c) {
            int g = word[c], k = g < 0 ? -g : g;
            if (pos == k)          { stepC[steps] = c; stepOver[steps] = g > 0; stepArcs[steps++] = arcs; arcs = 0; pos = k + 1; }
            else if (pos == k + 1) { stepC[steps] = c; stepOver[steps] = g < 0; stepArcs[steps++] = arcs; arcs = 0; pos = k; }
        }
        ++strands;
        if (pos == 1) break;
        ++arcs;
        x = pos;
    }
    if (strands != m) return false;
    int over[kMaxN], under[kMaxN];
    for (int idx = 0; idx < steps; ++idx) (stepOver[idx] ? over : under)[stepC[idx]] = idx + 1;
    K.n = n;
    for (int c = 0; c < n; ++c) K.x[c] = {word[c] > 0 ? 1 : -1, over[c], under[c]};
    for (int idx = 0; idx < steps; ++idx) K.rot[idx] = -stepArcs[idx];
    K.rot[steps] = 0;
    return true;
}

// ============================================================================
// 4. Gauss-Jordan mod q: (det, A^-1). False if A is singular at this point.
// ============================================================================
using Mat = u32[kMaxM][kMaxM];

bool detInv(Mat A, Mat Inv, int m, u32& det) {
    for (int r = 0; r < m; ++r) {
        for (int c = 0; c < m; ++c) Inv[r][c] = (r == c);
    }
    det = 1;
    for (int k = 0; k < m; ++k) {
        int p = k;
        while (p < m && A[p][k] == 0) ++p;
        if (p == m) return false;
        if (p != k) {
            for (int c = 0; c < m; ++c) { std::swap(A[k][c], A[p][c]); std::swap(Inv[k][c], Inv[p][c]); }
            det = neg(det);
        }
        u32 pv = A[k][k];
        det = mul(det, pv);
        u32 ip = inv(pv);
        for (int c = k; c < m; ++c) A[k][c] = mul(A[k][c], ip);
        for (int c = 0; c < m; ++c) Inv[k][c] = mul(Inv[k][c], ip);
        for (int r = 0; r < m; ++r) {
            if (r == k || A[r][k] == 0) continue;
            u32 f = A[r][k];
            for (int c = k; c < m; ++c) A[r][c] = sub(A[r][c], mul(f, A[k][c]));
            for (int c = 0; c < m; ++c) Inv[r][c] = sub(Inv[r][c], mul(f, Inv[k][c]));
        }
    }
    return true;
}

// ============================================================================
// 5. theta(D)(T1, T2) in F_q. False if T2 = 1 or T is a root of Delta (singular matrix).
// ============================================================================
bool thetaEval(const Knot& d, u32 T1, u32 T2, u32& theta, u32* deltas = nullptr) {
    const int n = d.n, m = 2 * n + 1;
    if (T2 == 1 || T1 == 0 || T2 == 0) return false;
    const u32 T3 = mul(T1, T2);
    const u32 half = static_cast<u32>((Q + 1) / 2);
    const u32 Ts[3] = {T1, T2, T3};

    Mat G[3];
    u32 Del = 1;
    const int ed = d.expDelta();
    for (int nu = 0; nu < 3; ++nu) {
        const u32 T = Ts[nu], Ti = inv(T);
        Mat A;
        for (int r = 0; r < m; ++r) for (int c = 0; c < m; ++c) A[r][c] = (r == c);
        for (int c = 0; c < n; ++c) {
            const u32 Tc = d.x[c].s > 0 ? T : Ti;
            const int i = d.x[c].i - 1, j = d.x[c].j - 1;
            A[i][i + 1] = sub(A[i][i + 1], Tc);
            A[i][j + 1] = sub(add(A[i][j + 1], Tc), 1);
            A[j][j + 1] = sub(A[j][j + 1], 1);
        }
        u32 det;
        if (!detInv(A, G[nu], m, det)) return false;
        const u32 dk = mul(ed >= 0 ? power(T, ed) : power(Ti, -ed), det);      // Delta(T) in F_q
        if (deltas) deltas[nu] = dk;
        Del = mul(Del, dk);
    }
    const Mat &G1 = G[0], &G2 = G[1], &G3 = G[2];

    // T^s per crossing and 1/(T2^s - 1)
    u32 t1s[kMaxN], t2s[kMaxN], t3s[kMaxN], it2m1[kMaxN];
    const u32 iT1 = inv(T1), iT2 = inv(T2), iT3 = inv(T3);
    const u32 iPos = inv(sub(T2, 1)), iNeg = inv(sub(iT2, 1));
    for (int c = 0; c < n; ++c) {
        const bool p = d.x[c].s > 0;
        t1s[c] = p ? T1 : iT1; t2s[c] = p ? T2 : iT2; t3s[c] = p ? T3 : iT3;
        it2m1[c] = p ? iPos : iNeg;
    }

    u32 th0 = 0;
    // ---- F1 (sum over crossings) ----
    for (int c = 0; c < n; ++c) {
        const int i = d.x[c].i - 1, j = d.x[c].j - 1;
        const u32 g1ii = G1[i][i], g2ii = G2[i][i], g3ii = G3[i][i];
        const u32 g1ji = G1[j][i], g2ji = G2[j][i], g3ji = G3[j][i];
        const u32 g2ij = G2[i][j], g2jj = G2[j][j], g3jj = G3[j][j];
        const u32 a = t1s[c], b = t2s[c], e = t3s[c];
        const u32 b1 = sub(b, 1), e1 = sub(e, 1), a1 = sub(a, 1);
        u32 poly = sub(half, g3ii);
        poly = add(poly, mul(mul(b, g1ii), g2ji));
        poly = sub(poly, mul(mul(b, g3jj), g2ji));
        poly = sub(poly, mul(mul(b1, g3ii), g2ji));
        poly = add(poly, mul(mul(e1, g2ji), g3ji));
        poly = sub(poly, mul(g1ii, g2jj));
        poly = add(poly, mul(add(g3ii, g3ii), g2jj));
        poly = add(poly, mul(g1ii, g3jj));
        poly = sub(poly, mul(g2ii, g3jj));
        u32 f1 = mul(mul(a1, b), add(sub(mul(g3jj, g1ji), mul(g2jj, g1ji)), mul(mul(b, g1ji), g2ji)));
        u32 inner = sub(1, mul(b, g1ii));
        inner = add(inner, g2ij);
        inner = add(inner, mul(sub(b, 2), g2jj));
        inner = sub(inner, mul(mul(a1, add(b, 1)), g1ji));
        u32 frac = add(f1, mul(mul(e1, g3ji), inner));
        u32 term = add(poly, mul(frac, it2m1[c]));
        th0 = d.x[c].s > 0 ? add(th0, term) : sub(th0, term);
    }
    // ---- F2 (sum over pairs of crossings) ----
    for (int c0 = 0; c0 < n; ++c0) {
        const int i0 = d.x[c0].i - 1, j0 = d.x[c0].j - 1;
        const u32 a0 = sub(t1s[c0], 1), tt = t2s[c0];
        for (int c1 = 0; c1 < n; ++c1) {
            const int i1 = d.x[c1].i - 1, j1 = d.x[c1].j - 1;
            u32 pref = mul(mul(a0, mul(sub(t3s[c1], 1), it2m1[c1])), mul(G1[j1][i0], G3[j0][i1]));
            if (pref == 0) continue;
            if (d.x[c1].s < 0) pref = neg(pref);
            u32 in = add(mul(tt, G2[i1][i0]), G2[j1][j0]);
            in = sub(sub(in, mul(tt, G2[j1][i0])), G2[i1][j0]);
            th0 = add(th0, mul(pref, in));
        }
    }
    // ---- F3 (sum over edges) ----
    for (int k = 0; k < m; ++k) th0 = add(th0, mul(fromInt(d.rot[k]), sub(G3[k][k], half)));

    theta = mul(th0, Del);
    return true;
}

// ----------------------------------------------------------------------------
// Evaluation points (public) and multipoint digest
// ----------------------------------------------------------------------------
constexpr int kMaxK = 16;
constexpr char kDomPoint[] = "holocrypto-v2/point";
constexpr char kDomTheta[] = "holocrypto-v2/theta8";
constexpr char kDomTheta3[] = "holocrypto-v3/alexander+theta";   // optional v3 digest (also binds Delta)

inline void putBE(uint8_t* p, u64 v) { for (int b = 0; b < 8; ++b) p[b] = static_cast<uint8_t>(v >> (56 - 8 * b)); }
inline u64 getBE(const uint8_t* p) { u64 v = 0; for (int b = 0; b < 8; ++b) v = (v << 8) | p[b]; return v; }

struct Params {
    int K = 8;
    u64 seed = 0x484F4C4F43525950ULL;          // "HOLOCRYP"
    u32 T1[kMaxK] = {}, T2[kMaxK] = {};
    void derive() {
        constexpr size_t dl = sizeof(kDomPoint) - 1;
        for (int k = 0; k < K; ++k) {
            uint8_t msg[dl + 9], h[32];
            std::memcpy(msg, kDomPoint, dl);
            putBE(msg + dl, seed);
            msg[dl + 8] = static_cast<uint8_t>(k);
            Sha256::hash(msg, sizeof(msg), h);
            T1[k] = static_cast<u32>(2 + getBE(h) % (Q - 3));
            T2[k] = static_cast<u32>(2 + getBE(h + 8) % (Q - 3));
        }
    }
};

// theta at the K points. False if the matrix is singular at any of them (the knot is rejected).
bool thetaMulti(const Knot& d, const Params& P, u32 th[kMaxK]) {
    for (int k = 0; k < P.K; ++k)
        if (!thetaEval(d, P.T1[k], P.T2[k], th[k])) return false;
    return true;
}

// SHA-256( "holocrypto-v2/theta8" || q || K || seed || theta_0 .. theta_{K-1} ), 8-byte big-endian integers
void holocryptoHash(const Params& P, const u32 th[kMaxK], uint8_t out[32]) {
    constexpr size_t dl = sizeof(kDomTheta) - 1;
    uint8_t msg[dl + 24 + 8 * kMaxK];
    std::memcpy(msg, kDomTheta, dl);
    putBE(msg + dl, Q);
    putBE(msg + dl + 8, static_cast<u64>(P.K));
    putBE(msg + dl + 16, P.seed);
    for (int k = 0; k < P.K; ++k) putBE(msg + dl + 24 + 8 * k, th[k]);
    Sha256::hash(msg, dl + 24 + 8 * P.K, out);
}

// v3 (optional, --verify --v3 only): SHA-256( "holocrypto-v3/alexander+theta" || q || K || seed ||
//   for each point k: Delta(T1_k), Delta(T2_k), Delta(T3_k), theta_k ), 8-byte big-endian integers.
// Same thetas as v2; the Alexander values are hashed too, so K # K* (theta = 0) no longer collapses onto the
// unknot digest. It does NOT separate knots that share Delta (e.g. Delta = 1: Conway, Kinoshita-Terasaka).
bool thetaMultiV3(const Knot& d, const Params& P, u32 th[kMaxK], u32 dl[kMaxK][3]) {
    for (int k = 0; k < P.K; ++k)
        if (!thetaEval(d, P.T1[k], P.T2[k], th[k], dl[k])) return false;
    return true;
}

void holocryptoHashV3(const Params& P, const u32 th[kMaxK], const u32 dl[kMaxK][3], uint8_t out[32]) {
    constexpr size_t dn = sizeof(kDomTheta3) - 1;
    uint8_t msg[dn + 24 + 32 * kMaxK];
    std::memcpy(msg, kDomTheta3, dn);
    putBE(msg + dn, Q);
    putBE(msg + dn + 8, static_cast<u64>(P.K));
    putBE(msg + dn + 16, P.seed);
    for (int k = 0; k < P.K; ++k) {
        uint8_t* w = msg + dn + 24 + 32 * k;
        putBE(w, dl[k][0]); putBE(w + 8, dl[k][1]); putBE(w + 16, dl[k][2]); putBE(w + 24, th[k]);
    }
    Sha256::hash(msg, dn + 24 + 32 * P.K, out);
}

// ============================================================================
// 6. Benchmark: all words of length L over {1,2,3,-1,-2,-3} (4 strands)
// ============================================================================
constexpr int kStrands = 4;
constexpr int8_t kAlphabet[6] = {1, 2, 3, -1, -2, -3};

struct Stats {
    u64 words = 0, links = 0, singular = 0, zeroAll = 0, zeroSome = 0, evals = 0, checksum = 0;
    void merge(const Stats& o) {
        words += o.words; links += o.links; singular += o.singular; zeroAll += o.zeroAll; zeroSome += o.zeroSome;
        evals += o.evals; checksum += o.checksum;
    }
};

void worker(int len, u64 lo, u64 hi, const Params& P, Stats& st) {
    int8_t w[16];
    Knot K;
    for (u64 idx = lo; idx < hi; ++idx) {
        u64 r = idx;
        for (int p = 0; p < len; ++p) { w[p] = kAlphabet[r % 6]; r /= 6; }
        ++st.words;
        if (!braidToKnot(w, len, kStrands, K)) { ++st.links; continue; }
        u32 th[kMaxK];
        if (!thetaMulti(K, P, th)) { ++st.singular; continue; }
        uint8_t h[32];
        holocryptoHash(P, th, h);
        u64 first; std::memcpy(&first, h, 8);
        st.checksum += first;                      // keeps the compiler from eliminating the computation
        int zeros = 0;
        for (int k = 0; k < P.K; ++k) zeros += (th[k] == 0);
        st.zeroAll += (zeros == P.K);              // theta = 0 at all points: key rejected by the filter
        st.zeroSome += (zeros > 0 && zeros < P.K); // zero at only some points (bad luck with the point)
        ++st.evals;
    }
}

Stats runLength(int len, int threads, const Params& P, double& secs) {
    u64 total = 1;
    for (int i = 0; i < len; ++i) total *= 6;
    std::vector<Stats> parts(threads);
    std::vector<std::thread> pool;
    auto t0 = std::chrono::steady_clock::now();
    for (int t = 0; t < threads; ++t) {
        u64 lo = total * t / threads, hi = total * (t + 1) / threads;
        pool.emplace_back(worker, len, lo, hi, std::cref(P), std::ref(parts[t]));
    }
    for (auto& th : pool) th.join();
    secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    Stats s;
    for (auto& p : parts) s.merge(p);
    return s;
}

// Strict integer list parser: every token must be a complete decimal integer within [lo, hi]; no silent
// wrap-around (atoi + int8_t cast turned "257" into 1 and "-256" into 0). Returns false and prints why.
bool parseIntList(const std::string& s, char sep, long lo, long hi, const char* what, std::vector<int>& out) {
    out.clear();
    std::string tok;
    std::istringstream is(s);
    while (std::getline(is, tok, sep)) {
        if (tok.empty()) continue;
        char* end = nullptr;
        errno = 0;
        long x = std::strtol(tok.c_str(), &end, 10);
        if (end == tok.c_str() || *end != '\0' || errno == ERANGE || x < lo || x > hi) {
            std::fprintf(stderr, "error: invalid %s '%s' (expected an integer in [%ld, %ld])\n", what, tok.c_str(), lo, hi);
            return false;
        }
        out.push_back(static_cast<int>(x));
    }
    return true;
}

int cmdVerify(const std::string& wordStr, const Params& P, bool v3) {
    std::vector<int> v;
    // letters must be sigma_k^{+-1} with 1 <= k <= kStrands-1 (0 and |k| >= kStrands are not generators)
    if (!parseIntList(wordStr, ' ', -(kStrands - 1), kStrands - 1, "braid letter", v)) return 2;
    for (int x : v) if (x == 0) { std::fprintf(stderr, "error: braid letter 0 is not a generator\n"); return 2; }
    if (v.size() > static_cast<size_t>(kMaxN)) {
        std::fprintf(stderr, "error: word has %zu letters; this build supports at most %d crossings\n", v.size(), kMaxN);
        return 2;
    }
    std::vector<int8_t> w(v.begin(), v.end());
    Knot K;
    if (!braidToKnot(w.data(), static_cast<int>(w.size()), kStrands, K)) { std::puts("link"); return 1; }
    u32 th[kMaxK], dl[kMaxK][3];
    if (!thetaMultiV3(K, P, th, dl)) { std::puts("singular"); return 1; }
    uint8_t h[32];
    if (v3) holocryptoHashV3(P, th, dl, h); else holocryptoHash(P, th, h);
    for (int k = 0; k < P.K; ++k) {
        std::printf("theta[%d]=%u  (T1=%u T2=%u)\n", k, th[k], P.T1[k], P.T2[k]);
        if (v3) std::printf("delta[%d]=%u %u %u\n", k, dl[k][0], dl[k][1], dl[k][2]);
    }
    std::printf("hash=%s\n", hex(h, 32).c_str());
    return 0;
}

int cmdSelftest() {
    struct { const char* in; const char* out; } tv[] = {
        {"", "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"},
        {"abc", "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"},
        {"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq",
         "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"}};
    bool ok = true;
    for (auto& t : tv) {
        uint8_t h[32];
        Sha256::hash(reinterpret_cast<const uint8_t*>(t.in), std::strlen(t.in), h);
        bool good = hex(h, 32) == t.out;
        ok &= good;
        std::printf("SHA-256(\"%.12s%s\") %s\n", t.in, std::strlen(t.in) > 12 ? "..." : "", good ? "OK" : "FAIL");
    }
    return ok ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    if (threads < 1) threads = 1;
    Params P;
    std::vector<int> lens = {5, 6, 7, 8};
    std::string verifyWord;
    bool v3 = false;                          // --v3: optional v3 digest, only meaningful with --verify
    bool haveVerify = false;                  // --verify "" must be an error, not a silent benchmark run
    for (int a = 1; a < argc; ++a) {
        std::string k = argv[a];
        auto next = [&]() -> const char* { return a + 1 < argc ? argv[++a] : ""; };
        std::vector<int> one;
        if (k == "-t") { if (!parseIntList(next(), ',', 1, 256, "thread count", one) || one.size() != 1) return 2; threads = one[0]; }
        else if (k == "-k") { if (!parseIntList(next(), ',', 1, kMaxK, "K", one) || one.size() != 1) return 2; P.K = one[0]; }
        else if (k == "--seed") {
            const char* t = next();
            char* end = nullptr;
            errno = 0;
            u64 sv = std::strtoull(t, &end, 0);
            // reject "", text, trailing junk, overflow and negative numbers (strtoull would wrap "-1" to 2^64-1)
            if (end == t || *end != '\0' || errno == ERANGE || std::strchr(t, '-')) {
                std::fprintf(stderr, "error: invalid seed '%s' (expected an unsigned 64-bit integer)\n", t);
                return 2;
            }
            P.seed = sv;
        }
        else if (k == "--lens") { if (!parseIntList(next(), ',', 1, kMaxN, "length", lens) || lens.empty()) return 2; }
        else if (k == "--verify") { verifyWord = next(); haveVerify = true; }
        else if (k == "--v3") v3 = true;
        else if (k == "--selftest") return cmdSelftest();
        else { std::fprintf(stderr, "error: unknown option '%s'\n", k.c_str()); return 2; }
    }
    P.derive();
    if (v3 && !haveVerify) { std::fprintf(stderr, "error: --v3 only applies to --verify\n"); return 2; }
    if (haveVerify) return cmdVerify(verifyWord, P, v3);

    std::printf("holocrypto_bench v2 | q=2^31-1 | K=%d points | seed=0x%llx | strands=%d | threads=%d\n", P.K,
                (unsigned long long)P.seed, kStrands, threads);
    std::printf("alphabet {1,2,3,-1,-2,-3}; iteration = closure -> theta at K points -> SHA-256\n\n");
    std::printf("%-4s %10s %9s %9s %9s %10s %12s %14s\n", "len", "words", "links", "singular", "knots", "time(s)",
                "words/s", "iterations/s");
    Stats all;
    double allSecs = 0;
    for (int len : lens) {
        double secs;
        Stats s = runLength(len, threads, P, secs);
        std::printf("%-4d %10llu %9llu %9llu %9llu %10.3f %12.0f %14.0f\n", len, (unsigned long long)s.words,
                    (unsigned long long)s.links, (unsigned long long)s.singular, (unsigned long long)s.evals, secs,
                    s.words / secs, secs > 0 ? s.evals / secs : 0.0);
        all.merge(s);
        allSecs += secs;
    }
    std::printf("\nTOTAL words=%llu knots evaluated=%llu (theta=0 at all %d points: %llu | at some: %llu) in %.3f s\n",
                (unsigned long long)all.words, (unsigned long long)all.evals, P.K, (unsigned long long)all.zeroAll,
                (unsigned long long)all.zeroSome, allSecs);
    std::printf("      -> %.0f iterations/s (%.2f us/iter amortised) | %.0f words/s\n", all.evals / allSecs,
                1e6 * allSecs / all.evals, all.words / allSecs);
    std::printf("checksum=%016llx\n", (unsigned long long)all.checksum);
    return 0;
}
