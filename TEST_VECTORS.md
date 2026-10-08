# Test vectors

Reference inputs and outputs for HOLO-CRYPTO (v2, multi-point). Everything below was produced by `holocrypto_bench.exe --verify` built from the sources in this repository (`zig c++ -O3 -march=native -std=c++17`), and **each vector was also recomputed with the independent Python reference** (`holocrypto_multipoint_ref.py`); all of them agree bit for bit.

Parameters unless stated otherwise: `q = 2147483647` (2^31 - 1), `K = 8`, seed `0x484F4C4F43525950`, 4 strands. Words are braid words over `{1, 2, 3, -1, -2, -3}` (`k` = sigma_k, `-k` = its inverse). Lines show `theta[k]` and the evaluation point `(T1, T2)` for each `k`, then the SHA-256 digest.

To reproduce the Conway and Kinoshita-Terasaka vectors, build the `kMaxN = 16` variant described in the README (section 5) and run `holocrypto_bench_n16.exe` instead.

## 1. Basic knot

5-7 crossing toy knot used throughout the README.

```text
$ holocrypto_bench.exe --verify "1 2 3 -1 2 1 3"
theta[0]=1883451254  (T1=1496904601 T2=1786645906)
theta[1]=260077083  (T1=735303977 T2=540787968)
theta[2]=1120809409  (T1=1585694247 T2=1330651864)
theta[3]=331051820  (T1=344924129 T2=330897688)
theta[4]=609041965  (T1=485100673 T2=294071181)
theta[5]=642464363  (T1=784128092 T2=174060674)
theta[6]=1355435651  (T1=2057504001 T2=2008740237)
theta[7]=1742479061  (T1=1052497489 T2=1421887965)
hash=b3afccc9aee25cfa5249a0089ac1341c9f843e4d2bdb945abc7c911c0549df31
```

## 2. Mirror image

All braid letters negated. Every theta becomes q - theta (q = 2147483647); the digest is unrelated.

```text
$ holocrypto_bench.exe --verify "-1 -2 -3 1 -2 -1 -3"
theta[0]=264032393  (T1=1496904601 T2=1786645906)
theta[1]=1887406564  (T1=735303977 T2=540787968)
theta[2]=1026674238  (T1=1585694247 T2=1330651864)
theta[3]=1816431827  (T1=344924129 T2=330897688)
theta[4]=1538441682  (T1=485100673 T2=294071181)
theta[5]=1505019284  (T1=784128092 T2=174060674)
theta[6]=792047996  (T1=2057504001 T2=2008740237)
theta[7]=405004586  (T1=1052497489 T2=1421887965)
hash=af3223d5dd5b165d371714a0ba7a52ca78cf0cc28e0eec2e9941dcb2704cbfa6
```

## 3. Five-crossing knot

Shortest legal length with a nonzero theta.

```text
$ holocrypto_bench.exe --verify "1 2 3 1 2"
theta[0]=1061173176  (T1=1496904601 T2=1786645906)
theta[1]=1209321409  (T1=735303977 T2=540787968)
theta[2]=544552224  (T1=1585694247 T2=1330651864)
theta[3]=881224373  (T1=344924129 T2=330897688)
theta[4]=34431997  (T1=485100673 T2=294071181)
theta[5]=520406002  (T1=784128092 T2=174060674)
theta[6]=478103607  (T1=2057504001 T2=2008740237)
theta[7]=1731653012  (T1=1052497489 T2=1421887965)
hash=17761f2b8ebe39ec012fa019cba7a93fce45ef85b171691cc849791dfd17a6ad
```

## 4. Its mirror

```text
$ holocrypto_bench.exe --verify "-1 -2 -3 -1 -2"
theta[0]=1086310471  (T1=1496904601 T2=1786645906)
theta[1]=938162238  (T1=735303977 T2=540787968)
theta[2]=1602931423  (T1=1585694247 T2=1330651864)
theta[3]=1266259274  (T1=344924129 T2=330897688)
theta[4]=2113051650  (T1=485100673 T2=294071181)
theta[5]=1627077645  (T1=784128092 T2=174060674)
theta[6]=1669380040  (T1=2057504001 T2=2008740237)
theta[7]=415830635  (T1=1052497489 T2=1421887965)
hash=05d50e2371f4d12fe62d1e1a0cb1b20ab9b583842ffbc3ef50bc9d60d038b06a
```

## 5. Equivalent diagram: sigma sigma^-1 inserted

Reidemeister II move on the basic knot. Must give the same digest as the first vector.

```text
$ holocrypto_bench.exe --verify "2 -2 1 2 3 -1 2 1 3"
theta[0]=1883451254  (T1=1496904601 T2=1786645906)
theta[1]=260077083  (T1=735303977 T2=540787968)
theta[2]=1120809409  (T1=1585694247 T2=1330651864)
theta[3]=331051820  (T1=344924129 T2=330897688)
theta[4]=609041965  (T1=485100673 T2=294071181)
theta[5]=642464363  (T1=784128092 T2=174060674)
theta[6]=1355435651  (T1=2057504001 T2=2008740237)
theta[7]=1742479061  (T1=1052497489 T2=1421887965)
hash=b3afccc9aee25cfa5249a0089ac1341c9f843e4d2bdb945abc7c911c0549df31
```

## 6. Equivalent diagram: cyclic rotation

Conjugation (cyclic rotation of the word). Must give the same digest as the first vector.

```text
$ holocrypto_bench.exe --verify "2 3 -1 2 1 3 1"
theta[0]=1883451254  (T1=1496904601 T2=1786645906)
theta[1]=260077083  (T1=735303977 T2=540787968)
theta[2]=1120809409  (T1=1585694247 T2=1330651864)
theta[3]=331051820  (T1=344924129 T2=330897688)
theta[4]=609041965  (T1=485100673 T2=294071181)
theta[5]=642464363  (T1=784128092 T2=174060674)
theta[6]=1355435651  (T1=2057504001 T2=2008740237)
theta[7]=1742479061  (T1=1052497489 T2=1421887965)
hash=b3afccc9aee25cfa5249a0089ac1341c9f843e4d2bdb945abc7c911c0549df31
```

## 7. Link (rejected)

Even length on 4 strands: the closure has several components, so the generator rejects it. No theta is computed.

```text
$ holocrypto_bench.exe --verify "1 2 3 -1 2 1"
link
```

## 8. Unknot (theta = 0)

Shortest legal length (n = m - 1) closes to the unknot. All 8 thetas are 0, so a key generator rejects it.

```text
$ holocrypto_bench.exe --verify "1 2 3"
theta[0]=0  (T1=1496904601 T2=1786645906)
theta[1]=0  (T1=735303977 T2=540787968)
theta[2]=0  (T1=1585694247 T2=1330651864)
theta[3]=0  (T1=344924129 T2=330897688)
theta[4]=0  (T1=485100673 T2=294071181)
theta[5]=0  (T1=784128092 T2=174060674)
theta[6]=0  (T1=2057504001 T2=2008740237)
theta[7]=0  (T1=1052497489 T2=1421887965)
hash=0a61a9a29d9315849958810aebc49c052249fa9946f88f4ccc2774e5902a96d4
```

## 9. Amphichiral / trivial (theta = 0)

Structural zero at all 8 points (the README's section 3.5 example).

```text
$ holocrypto_bench.exe --verify "1 -2 3 1 -2 3 -1"
theta[0]=0  (T1=1496904601 T2=1786645906)
theta[1]=0  (T1=735303977 T2=540787968)
theta[2]=0  (T1=1585694247 T2=1330651864)
theta[3]=0  (T1=344924129 T2=330897688)
theta[4]=0  (T1=485100673 T2=294071181)
theta[5]=0  (T1=784128092 T2=174060674)
theta[6]=0  (T1=2057504001 T2=2008740237)
theta[7]=0  (T1=1052497489 T2=1421887965)
hash=0a61a9a29d9315849958810aebc49c052249fa9946f88f4ccc2774e5902a96d4
```

## 10. Conway knot 11n34

11 crossings. Same Alexander polynomial and volume as Kinoshita-Terasaka; the two are compared on the kMaxN = 16 build.

```text
$ holocrypto_bench_n16.exe --verify "2 2 2 1 -3 -2 -2 1 -2 1 -3"
theta[0]=2099079240  (T1=1496904601 T2=1786645906)
theta[1]=926803104  (T1=735303977 T2=540787968)
theta[2]=1921730829  (T1=1585694247 T2=1330651864)
theta[3]=1405934459  (T1=344924129 T2=330897688)
theta[4]=1253585632  (T1=485100673 T2=294071181)
theta[5]=206900144  (T1=784128092 T2=174060674)
theta[6]=135011621  (T1=2057504001 T2=2008740237)
theta[7]=1311871219  (T1=1052497489 T2=1421887965)
hash=3e09f21ecc6f12dfa46a2647d714e3f668431bf923a4edc26474f7ba59e4bbb6
```

## 11. Kinoshita-Terasaka knot 11n42

13 crossings: needs the kMaxN = 16 build. Same Alexander polynomial and volume as Conway; all 8 thetas differ.

```text
$ holocrypto_bench_n16.exe --verify "1 1 1 3 3 2 -3 -1 -1 2 -1 -3 -2"
theta[0]=1764820622  (T1=1496904601 T2=1786645906)
theta[1]=1767663427  (T1=735303977 T2=540787968)
theta[2]=465091019  (T1=1585694247 T2=1330651864)
theta[3]=830736253  (T1=344924129 T2=330897688)
theta[4]=1008389583  (T1=485100673 T2=294071181)
theta[5]=63131087  (T1=784128092 T2=174060674)
theta[6]=2108633172  (T1=2057504001 T2=2008740237)
theta[7]=145495178  (T1=1052497489 T2=1421887965)
hash=5289aa05b81d57c5015f77af6559c6de4d590d118bac601df9964eda98aac608
```

## 12. Different public seed

Basic knot, --seed 0xC0FFEE: different 8 points, different thetas and digest.

```text
$ holocrypto_bench.exe --verify "1 2 3 -1 2 1 3" --seed 0xC0FFEE
theta[0]=1452855192  (T1=905818359 T2=490686565)
theta[1]=1267844968  (T1=1749125137 T2=412755550)
theta[2]=1706380425  (T1=592249859 T2=2063778292)
theta[3]=1111221291  (T1=1805661527 T2=1915977015)
theta[4]=572703196  (T1=1753423131 T2=1371586812)
theta[5]=122966231  (T1=348699552 T2=1970932708)
theta[6]=171674692  (T1=802336147 T2=204602502)
theta[7]=1021932365  (T1=145612200 T2=1380473798)
hash=66ddee70ecac08f4441225c30fcf6df8bb984eeb0f1099fa90268e37a6b50fc8
```

## 13. Single point (K = 1)

Basic knot, -k 1 (the v1-style single-point evaluation, but with the v2 domain tags).

```text
$ holocrypto_bench.exe --verify "1 2 3 -1 2 1 3" -k 1
theta[0]=1883451254  (T1=1496904601 T2=1786645906)
hash=b79fbd4713bd9d8a6c7b5a6c0e2fd1daf3f8c7faa6defc9ffb3e27d1519d5c62
```

## 14. K = 4

```text
$ holocrypto_bench.exe --verify "1 2 3 -1 2 1 3" -k 4
theta[0]=1883451254  (T1=1496904601 T2=1786645906)
theta[1]=260077083  (T1=735303977 T2=540787968)
theta[2]=1120809409  (T1=1585694247 T2=1330651864)
theta[3]=331051820  (T1=344924129 T2=330897688)
hash=1572d81bb91191565f8fac0d5047711c980c0e8ce7f27d2ec99d946bf7c9bde7
```

## 15. K = 16

Maximum K supported by the build.

```text
$ holocrypto_bench.exe --verify "1 2 3 -1 2 1 3" -k 16
theta[0]=1883451254  (T1=1496904601 T2=1786645906)
theta[1]=260077083  (T1=735303977 T2=540787968)
theta[2]=1120809409  (T1=1585694247 T2=1330651864)
theta[3]=331051820  (T1=344924129 T2=330897688)
theta[4]=609041965  (T1=485100673 T2=294071181)
theta[5]=642464363  (T1=784128092 T2=174060674)
theta[6]=1355435651  (T1=2057504001 T2=2008740237)
theta[7]=1742479061  (T1=1052497489 T2=1421887965)
theta[8]=1109205108  (T1=334499100 T2=1887983737)
theta[9]=1813248934  (T1=1069239362 T2=268667058)
theta[10]=1692475624  (T1=177325026 T2=1850688044)
theta[11]=159205733  (T1=678944130 T2=1528148831)
theta[12]=43303691  (T1=1538718600 T2=400585448)
theta[13]=1133042524  (T1=11882424 T2=43942400)
theta[14]=1350742795  (T1=2062565434 T2=159181873)
theta[15]=1575970006  (T1=377114959 T2=1751274819)
hash=cf911bfca2310529f729da5c448eb43521ea84adbd57c1b199efac59fa0d09d7
```

## Consistency checks that hold for these vectors

- Vectors 1 and 2: `theta_mirror[k] = q - theta[k]` for every `k`.
- Vectors 1, 5 and 6 are three diagrams of the same knot and give the same digest.
- Vector 7 is rejected as a link; vectors 8 and 9 are rejected because `theta = 0` at all 8 points.
- Vectors 10 and 11 (Conway / Kinoshita-Terasaka) have identical Alexander polynomial (`1`) and hyperbolic volume, yet all 8 thetas differ (`python mutant_test.py` checks this and the invariance under 60 rewritings of each word).
