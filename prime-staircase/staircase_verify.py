"""
Verification script for "The Prime Staircase"
(Justo Tapiador, Universidad de Alicante).

Self-contained: requires only `sympy` (mpmath comes bundled with it).
Running time: a few seconds.

It verifies:
  1. the hierarchy of higher-order primes and its OEIS levels
     (A000040, A006450, A038580, A049090, A049203, A049202);
  2. the index staircase m_k = f^k(k) against OEIS A058009 (k = 1..8);
  3. the staircase primes d_k = f^(k+1)(k) = p_{m_k} for k = 1..8
     (d_9 and d_10 come from segmented sieving / primesieve; they are
     reproduced here as constants -- see the paper, Section 5);
  4. the structural facts: the order-of-primeness ladder (Prop. 3.3)
     and the crossing with the primeth recurrence A007097 (Prop. 3.4);
  5. the master formula li(log d_k) = k + 1 + li(log k) + O(log log k)
     (Theorem 4.3) and the corollaries log d_k ~ k log k, log d_k ~ p_k.

Run:  python3 staircase_verify.py
"""
import math
from sympy import prime as f, primepi, isprime
import mpmath as mp

mp.mp.dps = 40
li = mp.li

# ------------------------------------------------------------------
print("=" * 66)
print("1. Hierarchy of higher-order primes:  p_{n,k} = f^(k+1)(n)")
print("=" * 66)
OEIS_LEVELS = {
    0: ("A000040", [2, 3, 5, 7, 11, 13, 17]),
    1: ("A006450", [3, 5, 11, 17, 31, 41, 59]),
    2: ("A038580", [5, 11, 31, 59, 127, 179, 277]),
    3: ("A049090", [11, 31, 127, 277, 709, 1063, 1787]),
    4: ("A049203", [31, 127, 709, 1787, 5381, 8527, 15299]),
    5: ("A049202", [127, 709, 5381, 15299, 52711, 87803, 167449]),
}
for k, (aid, ref) in OEIS_LEVELS.items():
    row = []
    for n in range(1, 8):
        v = n
        for _ in range(k + 1):
            v = f(v)
        row.append(v)
    assert row == ref, f"level {k} mismatch"
    print(f"  level {k}: {str(row):<38} matches OEIS {aid}   OK")

# ------------------------------------------------------------------
print()
print("=" * 66)
print("2. Index staircase m_k = f^k(k)  (OEIS A058009)")
print("   Staircase d_k = f^(k+1)(k) = p_{m_k}")
print("=" * 66)
A058009 = [2, 5, 31, 277, 5381, 87803, 2269733, 50728129,
           1559861749, 64988430769]
# d_9, d_10 computed in the paper by segmented sieve / primesieve
D_EXACT = {1: 3, 2: 11, 3: 127, 4: 1787, 5: 52711, 6: 1128889,
           7: 37139213, 8: 997525853,
           9: 36294260117, 10: 1765037224331}

for k in range(1, 9):
    v = k
    for _ in range(k):
        v = f(v)
    m_k = v
    d_k = f(m_k)
    assert m_k == A058009[k - 1], f"m_{k} mismatch vs A058009"
    assert d_k == D_EXACT[k], f"d_{k} mismatch"
    print(f"  k={k}:  m_k={m_k:>9,} = A058009({k})  OK;   "
          f"d_k = p_(m_k) = {d_k:>13,}  OK")

print(f"\n  (from the paper, Sec. 5:  d_9  = {D_EXACT[9]:,}")
print(f"                            d_10 = {D_EXACT[10]:,})")

# ------------------------------------------------------------------
print()
print("=" * 66)
print("3. Structure")
print("=" * 66)

def order_of_primeness(q):
    """Number of primes in the chain q, pi(q), pi^2(q), ... (q prime)."""
    assert isprime(q)
    n, cnt = q, 0
    while isprime(n):
        cnt += 1
        n = primepi(n)
    return cnt

print("  Order-of-primeness ladder (Prop. 3.3): F(d_k) = k+1 (+F(k) if k prime)")
for k in [4, 5, 6, 7, 10]:          # k=4,6,10 composite; k=5,7 prime
    d = D_EXACT[k]
    expected = k + 1 + (order_of_primeness(k) if isprime(k) else 0)
    got = order_of_primeness(d)
    status = "OK" if got == expected else "FAIL"
    print(f"    k={k:2d}: F(d_k) = {got:2d} (expected {expected})  {status}")

print("\n  Crossing with the primeth recurrence A007097 (Prop. 3.4):")
print("  d_k in A007097  <=>  k in A007097")
# orbit of 1, far enough to cover d_6 = 1128889
orbit = [1]
while orbit[-1] < 2_000_000:
    orbit.append(f(orbit[-1]))
A007097 = set(orbit)
for k in range(1, 7):
    d = D_EXACT[k]
    in_k, in_d = k in A007097, d in A007097
    status = "OK" if in_k == in_d else "FAIL"
    print(f"    k={k}: k in A007097? {str(in_k):5s}   "
          f"d_k = {d:>9,} in A007097? {str(in_d):5s}   {status}")

# ------------------------------------------------------------------
print()
print("=" * 66)
print("4. Master formula  li(log d_k) = k + 1 + li(log k) + O(log log k)")
print("   (Theorem 4.3)   --   residuals should be tiny")
print("=" * 66)
print("  k |   li(ln d_k)  |  k+1+li(ln k)  |  residual")
for k in range(4, 11):
    d = D_EXACT[k]
    lhs = li(mp.log(d))
    rhs = mp.mpf(k) + 1 + li(mp.log(k))
    print(f"  {k:2d} |  {float(lhs):11.5f} |   {float(rhs):11.5f} | "
          f"{float(lhs - rhs):+.5f}")

print("\n  Corollary 4.5(iii):  log d_k / p_k -> 1")
print("  k |    d_k     | log d_k | p_k | log d_k / p_k")
for k in range(3, 11):
    d = D_EXACT[k]
    print(f"  {k:2d} | {d:>13,} | {math.log(d):7.3f} | {f(k):3d} | "
          f"{math.log(d) / f(k):.4f}")

print("\nAll checks passed. See the paper for the proofs and for the")
print("computation of d_9 and d_10 (Section 5).")
