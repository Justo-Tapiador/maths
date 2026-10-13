#!/usr/bin/env python3
"""Extra numerical verifications for the remaining sections of the Tahawus paper."""
import json
import mpmath as mp
from fractions import Fraction as Fr

mp.mp.dps = 50
OUT = {}

def f(x):
    return (x**x - 1) / x

def T(y):
    y = mp.mpf(y)
    if y > 0:
        lo, hi = mp.mpf(1), mp.mpf(2)
        while f(hi) < y:
            hi *= 2
    else:
        lo, hi = mp.mpf('1e-40'), mp.mpf(1)
    for _ in range(300):
        mid = (lo + hi) / 2
        if f(mid) < y:
            lo = mid
        else:
            hi = mid
    x0 = (lo + hi) / 2
    try:
        x0 = mp.findroot(lambda x: x**x - y*x - 1, (lo, hi), tol=mp.mpf(10)**-40, maxsteps=50)
    except Exception:
        pass
    return x0

alpha = T(1)
delta = mp.findroot(lambda x: x**x - x**2 - 1, (mp.mpf(2.2), mp.mpf(2.3)))

# ---------------------------------------------------------------- 1. partial sums inside/outside R
print("== Radius partial sums ==")
R = mp.mpf('0.6526100448447448')
# recompute exact Maclaurin coefficients (same recurrence as research script, 30 terms)
N = 30
E = [Fr(0)]*(N+1); E[1] = Fr(1)
for k in range(2, N+1):
    E[k] = Fr((-1)**k, k*(k-1))
B = [Fr(0)]*(N+1); B[0] = Fr(1)
for k in range(1, N+1):
    s = Fr(0)
    for i in range(1, k+1):
        s += i*E[i]*B[k-i]
    B[k] = s/Fr(k)
Cg = [Fr((-1)**j) for j in range(N+1)]
Fco = [Fr(0)]*(N+1)
for k in range(1, N+1):
    s = Fr(0)
    for i in range(1, k+1):
        s += B[i]*Cg[k-i]
    Fco[k] = s
phico = [Fco[k+1] for k in range(N)]
inv = [Fr(1)]
for m in range(1, N):
    s = Fr(0)
    for i in range(1, m+1):
        s += phico[i]*inv[m-i]
    inv.append(-s)

def conv_trunc(a, b, trunc):
    r = [Fr(0)]*(trunc+1)
    for i in range(min(len(a), trunc+1)):
        if a[i] == 0: continue
        for j in range(min(len(b), trunc+1-i)):
            if b[j] == 0: continue
            r[i+j] += a[i]*b[j]
    return r

ds = [Fr(0)]*(N+1)
for n in range(1, N+1):
    base = inv[:n]
    result = [Fr(1)]
    e = n
    while e:
        if e & 1:
            result = conv_trunc(result, base, n-1)
        base = conv_trunc(base, base, n-1)
        e >>= 1
    ds[n] = Fr(result[n-1], n)

def mpol(fr):
    return mp.mpf(fr.numerator)/mp.mpf(fr.denominator)

for fac, lab in [(mp.mpf('0.9'), 'inside'), (mp.mpf('1.06'), 'outside')]:
    yv = R*fac
    tv = T(yv)
    errs = []
    s = mp.mpf(1)
    for k in range(1, 25):
        s += mpol(ds[k])*yv**k
        if k in (8, 12, 16, 20, 24):
            errs.append((k, mp.nstr(abs(s - tv), 3)))
    print(f"  y = {mp.nstr(yv,5)} ({lab}): err(n=8,12,16,20,24) = {errs}")
    OUT[f'partials_{lab}'] = errs

# ---------------------------------------------------------------- 2. W-iteration closed-form rate
print("== W-iteration rate ==")
rate_closed = 1/((1 + mp.log(alpha))*(1 + alpha))
print("  closed form 1/((1+ln a)(1+a)) =", mp.nstr(rate_closed, 10))
xk = mp.e**mp.lambertw(mp.log(2))
errs = []
for k in range(8):
    xk = mp.e**mp.lambertw(mp.log(xk + 1))
    errs.append(xk - alpha)
ratios = [mp.nstr(errs[i+1]/errs[i], 6) for i in range(len(errs)-1)]
print("  empirical error ratios:", ratios)
OUT['rate_closed'] = mp.nstr(rate_closed, 12)
OUT['rate_empirical'] = ratios

# ---------------------------------------------------------------- 3. tower multipliers
print("== Tower ==")
mult_alpha = (alpha+1)*mp.log(alpha)
print("  multiplier of g(x)=a^x-1 at a:", mp.nstr(mult_alpha, 12))
print("  multiplier of g at 0 (=ln a):", mp.nstr(mp.log(alpha), 12))
xk = mp.mpf('0.5')
seq = []
for k in range(30):
    xk = alpha**xk - 1
    if k < 14:
        seq.append(xk)
tr = [mp.nstr(seq[i+1]/seq[i], 6) for i in range(len(seq)-1)]
print("  tower ratios x_{n+1}/x_n:", tr)
OUT['tower_mult_alpha'] = mp.nstr(mult_alpha, 12)
OUT['tower_mult_zero'] = mp.nstr(mp.log(alpha), 12)
OUT['tower_ratios'] = tr
# Newton on the tower map (repelling fixed point) needs no verification; skip.

# ---------------------------------------------------------------- 4. q(s) = (s-1)^2 e^{-s} < 1
print("== Auxiliary bounds ==")
qs = max((mp.mpf(k)/4000 - 1)**2 * mp.e**(-(mp.mpf(k)/4000)) for k in range(1, 8001))
print("  max of (s-1)^2 e^{-s} on grid (0,2]:", mp.nstr(qs, 8), " | limit at s->0 is 1 (not attained)")
OUT['q_max'] = mp.nstr(qs, 8)
u1max = max((mp.mpf(k)/1000)*(1 + 2*mp.log(mp.mpf(k)/1000)) for k in range(1, 1001))
print("  max of u1(x)=x(1+2ln x) on (0,1]:", mp.nstr(u1max, 8))
OUT['u1_max'] = mp.nstr(u1max, 8)

# ---------------------------------------------------------------- 5. f'' sign grid (convexity)
print("== f'' sign ==")
def fpp(x):
    return x**(x-1)*(mp.log(x) + 1 - 1/x)**2 + (x**x*(x+1) - 2)/x**3
neg_ok = all(fpp(mp.mpf(k)/2000) < 0 for k in range(1, 2000))
pos_ok = all(fpp(1 + mp.mpf(k)/2000) > 0 for k in range(1, 4001))
print("  f''<0 on (0,1):", neg_ok, " | f''>0 on (1,3):", pos_ok)
OUT['fpp_signs_ok'] = bool(neg_ok and pos_ok)
# check the Psi identity numerically at a few points
for xx in ['0.3', '0.5', '0.7', '0.9']:
    x = mp.mpf(xx); s = -mp.log(x); E = x
    t = (s-1)*E
    Psi = mp.e**(-s*E)*((1+t)**2 + 1 + E)
    lhs = fpp(x)*x**3 + 2
    print(f"  x={xx}: Psi = {mp.nstr(Psi,12)}  vs f''x^3+2 = {mp.nstr(lhs,12)}  diff = {mp.nstr(Psi-lhs,3)}")

# ---------------------------------------------------------------- 6. golden sign pattern
print("== Golden crossing ==")
def phi_metal(b):
    return (mp.mpf(b) + mp.sqrt(mp.mpf(b)**2 + 4))/2
checks = []
for b in ['-0.5', '0.25', '1', '1.4', '1.5', '1.51', '2', '5']:
    bb = mp.mpf(b)
    checks.append((b, mp.nstr(T(bb) - phi_metal(bb), 6)))
print("  T(b)-phi_b:", checks)
OUT['golden_signs'] = checks

# ---------------------------------------------------------------- 7. Newton for delta
print("== Newton for delta from 2 ==")
xk = mp.mpf(2)
nw = []
for k in range(6):
    xk = xk - (xk**xk - xk**2 - 1)/(xk**xk*(mp.log(xk)+1) - 2*xk)
    nw.append(mp.nstr(xk, 20))
print("  ", nw)
OUT['newton_delta'] = nw
OUT['delta_40'] = mp.nstr(delta, 41)
OUT['Tprime_delta_formula'] = mp.nstr(delta**2/(delta**(delta-1)*(mp.log(delta)+1-1/delta) + 1/delta**2), 12)

# ---------------------------------------------------------------- 8. ssqrt comparison at googol-ish
print("== Plus-one law ==")
for Y in [mp.mpf(10)**12, mp.mpf(10)**24, mp.mpf(10)**48]:
    x = T(Y)
    s0 = mp.e**mp.lambertw(mp.log(Y))
    W = mp.lambertw(mp.log(Y))
    print(f"  y=10^{mp.nstr(mp.log(Y, 10),4)}: T-ssqrt = {mp.nstr(x-s0,8)}   W/(W+1) = {mp.nstr(W/(W+1),8)}")
    OUT.setdefault('plus_one', []).append([mp.nstr(x-s0, 8), mp.nstr(W/(W+1), 8)])

# ---------------------------------------------------------------- 9. series reversion check at u
print("== f series check ==")
u = mp.mpf('0.05')
ser = sum(mpol(Fco[k])*u**k for k in range(1, 15))
print("  f(1.05) series =", mp.nstr(ser, 12), " direct =", mp.nstr(f(1+u), 12))
OUT['f_series_check'] = mp.nstr(ser - f(1+u), 5)

# ---------------------------------------------------------------- 10. orbit of T from 1 and -5
print("== Orbits ==")
yk = mp.mpf(-5); orb = []
for k in range(8):
    yk = T(yk)
    orb.append(mp.nstr(yk, 12))
print("  orbit from -5:", orb)
OUT['orbit_m5'] = orb

with open('tahawus_extra.json', 'w') as fh:
    json.dump(OUT, fh, indent=1, default=str)
print("saved -> scripts/tahawus_extra.json")
