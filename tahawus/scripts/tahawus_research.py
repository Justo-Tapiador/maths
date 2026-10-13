#!/usr/bin/env python3
"""
Computational research on the Tahawus function  T(y) = f^{-1}(y),  f(x) = (x^x - 1)/x.
Produces all numerical facts used in the paper; saves results to tahawus_results.json.
"""
import json, sys, math, cmath
import mpmath as mp
import sympy as sp

mp.mp.dps = 80
OUT = {}

def f(x):
    return (x**x - 1) / x

def T(y):
    """Robust inverse of f on (0,inf)."""
    y = mp.mpf(y)
    if y > 0:
        lo, hi = mp.mpf(1), mp.mpf(2)
        while f(hi) < y:
            hi *= 2
    else:
        # T(y) in (0,1]; for very negative y, bisect in t = ln x
        if y < -5:
            tlo, thi = y - 2, mp.mpf(0)
            g = lambda t: (mp.e**(t * mp.e**t) - 1) * mp.e**(-t) - y
            for _ in range(300):
                tm = (tlo + thi) / 2
                if g(tm) < 0:
                    tlo = tm
                else:
                    thi = tm
            return mp.e**((tlo + thi) / 2)
        lo, hi = mp.mpf('1e-40'), mp.mpf(1)
    for _ in range(400):
        mid = (lo + hi) / 2
        if f(mid) < y:
            lo = mid
        else:
            hi = mid
    x0 = (lo + hi) / 2
    try:
        x0 = mp.findroot(lambda x: x**x - y * x - 1, (lo, hi),
                         tol=mp.mpf(10)**-70, maxsteps=60)
    except Exception:
        pass
    return x0

def phi_metal(b):
    return (mp.mpf(b) + mp.sqrt(mp.mpf(b)**2 + 4)) / 2

# ---------------------------------------------------------------- constants
alpha = T(1)      # Tahawus constant  (root x^x = x+1)
beta  = T(-1)     # root x^x = 1-x
v     = 1/beta - 1     # A085846
c1    = 1/beta         # A169862 = v + 1
delta = mp.findroot(lambda x: x**x - x**2 - 1, (mp.mpf(2.2), mp.mpf(2.3)),
                    tol=mp.mpf(10)**-70, maxsteps=80)     # fixed point of T
T2 = T(2)

OUT['alpha'] = mp.nstr(alpha, 75)
OUT['beta'] = mp.nstr(beta, 75)
OUT['v'] = mp.nstr(v, 75)
OUT['c1'] = mp.nstr(c1, 75)
OUT['delta'] = mp.nstr(delta, 60)
OUT['T2'] = mp.nstr(T2, 60)

print("== Constants ==")
print("alpha = T(1)  =", OUT['alpha'])
print("beta  = T(-1) =", OUT['beta'])
print("v=1/beta-1    =", OUT['v'])
print("c1=1/beta     =", OUT['c1'])
print("delta (fixpt) =", OUT['delta'])
print("T(2)          =", OUT['T2'])

# --- verify against OEIS digit strings (A124930, A085846) ---
oeis_a124930 = ("1.7767750400970546974797307440387567486374110343292961390843"
                "74015273118658932824770702072786151313523630 09".replace(" ", ""))
oeis_a085846 = ("2.2931662874118610315080282912508058643722572903271212485377"
                "103961685064880091577436290420138048282 566 1".replace(" ", ""))
blog_alpha = "1.7767750400970546974797307440387567486374110343292961390843740"
blog_beta  = "0.3036591270299660512450180"
def digits_agree(mine, ref, nn):
    sa = mp.nstr(mine, nn+2).replace('.', '', 1)
    sr = ref.replace('.', '', 1).replace(' ', '')
    n2 = min(nn, len(sa), len(sr))
    return sa[:n2] == sr[:n2]
for name, mine, ref, nn in [("alpha_vs_OEIS", alpha, oeis_a124930, 70),
                            ("alpha_vs_blog", alpha, blog_alpha, 65),
                            ("v_vs_OEIS", v, oeis_a085846, 70),
                            ("beta_vs_blog", beta, blog_beta, 22)]:
    ok = digits_agree(mine, ref, nn)
    print(f"CHECK {name}: {'OK' if ok else 'MISMATCH'}")
    OUT['check_' + name] = bool(ok)

# --- defining-equation checks ---
OUT['alpha_check'] = mp.nstr(alpha**alpha - alpha - 1, 5)
OUT['beta_check']  = mp.nstr(beta**beta + beta - 1, 5)
OUT['v_check']     = mp.nstr(v**(1 + 1/v) - v - 1, 5)          # v^(1+1/v)=v+1
OUT['c_check']     = mp.nstr(c1**(c1-1) - (c1-1)**c1, 5)       # c^(c-1)=(c-1)^c
OUT['xth1_over_x_pair'] = mp.nstr(v**(1/v) - (v+1)**(1/(v+1)), 5)
print("alpha^alpha-alpha-1 =", OUT['alpha_check'], "| beta^beta+beta-1 =", OUT['beta_check'])
print("v^(1+1/v)-v-1       =", OUT['v_check'], "| c^(c-1)-(c-1)^c =", OUT['c_check'])

# A144184 family: root of x^(1+1/x) = x+2
w144 = mp.findroot(lambda x: x**(1 + 1/x) - x - 2, (mp.mpf(5), mp.mpf(6)),
                   tol=mp.mpf(10)**-60)
OUT['A144184'] = mp.nstr(w144, 55)
OUT['A144184_check'] = mp.nstr(w144**(1 + 1/w144) - w144 - 2, 5)
print("A144184 (x^(1+1/x)=x+2) =", OUT['A144184'])

# golden-ratio crossing point of x^x and x^(1+1/x):  phi
phi = (1 + mp.sqrt(5)) / 2
OUT['phi_cross_check'] = mp.nstr(phi**phi - phi**(1 + 1/phi), 5)
print("phi^phi - phi^(1+1/phi) =", OUT['phi_cross_check'])

# ---------------------------------------------------------------- Maclaurin series (sympy)
print("\n== Maclaurin series of T at 0 (Lagrange-Buermann, exact rationals) ==")
N = 60
from fractions import Fraction as Fr
# E(u) = (1+u) ln(1+u) = u + sum_{k>=2} (-1)^k u^k/(k(k-1))
E = [Fr(0)]*(N+1)
E[1] = Fr(1)
for k in range(2, N+1):
    E[k] = Fr((-1)**k, k*(k-1))
# x^x = exp(E):  B_0 = 1;  k B_k = sum_{i=1..k} i E_i B_{k-i}
B = [Fr(0)]*(N+1); B[0] = Fr(1)
for k in range(1, N+1):
    s = Fr(0)
    for i in range(1, k+1):
        s += i*E[i]*B[k-i]
    B[k] = s/Fr(k)
# f(1+u) = (x^x - 1)/(1+u) = (B - B_0) conv (1,-1,1,-1,...)
Cg = [Fr((-1)**j) for j in range(N+1)]
Fco = [Fr(0)]*(N+1)
for k in range(1, N+1):
    s = Fr(0)
    for i in range(1, k+1):
        s += B[i]*Cg[k-i]
    Fco[k] = s
print("f(1+u) = u +", [str(Fco[k]) for k in range(2, 13)], "...")
OUT['F_coeffs'] = [str(Fco[k]) for k in range(0, 13)]
assert Fco[1] == 1 and Fco[2] == 0 and Fco[3] == Fr(1,2) and Fco[4] == Fr(-1,6)

# phi(t) = F(t)/t = sum phi[k] t^k,  phi[0] = 1
phico = [Fco[k+1] for k in range(N)]
# series inverse inv = phi^{-1} (constant term 1): classic recurrence
inv = [Fr(1)]
for m in range(1, N):
    s = Fr(0)
    for i in range(1, m+1):
        s += phico[i]*inv[m-i]
    inv.append(-s)

def conv_trunc(a, b, trunc):
    """convolution of coefficient lists a,b truncated to degree trunc (inclusive)"""
    r = [Fr(0)]*(trunc+1)
    for i in range(min(len(a), trunc+1)):
        if a[i] == 0: continue
        for j in range(min(len(b), trunc+1-i)):
            if b[j] == 0: continue
            r[i+j] += a[i]*b[j]
    return r

# d_n = (1/n) [t^{n-1}] phi(t)^{-n}
ds = [Fr(0)]*(N+1)
for n in range(1, N+1):
    # psi^n via binary exponentiation, truncated at degree n-1
    base = inv[:n]              # only need up to degree n-1
    result = [Fr(1)]
    e = n
    while e:
        if e & 1:
            result = conv_trunc(result, base, n-1)
        base = conv_trunc(base, base, n-1)
        e >>= 1
    ds[n] = Fr(result[n-1], n)
OUT['maclaurin_exact'] = [str(ds[k]) for k in range(1, 13)]
OUT['maclaurin_dec'] = [float(ds[k]) for k in range(1, N+1)]
print("T(y) = 1 +", [str(ds[k]) for k in range(1, 13)], "...")

# verify numerically against real T at several points
for ytest in ['0.15', '0.3', '0.45']:
    yv = mp.mpf(ytest)
    series_val = mp.mpf(1) + sum(mp.mpf(ds[k].numerator)/mp.mpf(ds[k].denominator)*yv**k
                                for k in range(1, 13))
    OUT[f'maclaurin_check_{ytest}'] = mp.nstr(series_val - T(yv), 8)
    print(f"  series({ytest}) - T({ytest}) =", OUT[f'maclaurin_check_{ytest}'])

# radius of convergence estimate from coefficients
ratios = []
for n in range(10, N):
    if ds[n] != 0 and ds[n+1] != 0:
        ratios.append(abs(float(ds[n])/float(ds[n+1])))
OUT['radius_ratios'] = ratios
OUT['root_test'] = [abs(float(ds[n]))**(-1.0/n) for n in range(12, N) if ds[n] != 0]
print("ratio estimates |b_n/b_{n+1}| (n=10..59):", [round(r, 4) for r in ratios[:8]], "...")
print("root estimates |b_n|^{-1/n} (n=45..59):", [round(r, 5) for r in OUT['root_test'][-12:]])

# ---------------------------------------------------------------- nearest singularities
print("\n== Nearest singularities of T: analysis via Maclaurin coefficients ==")
def Fc(x):
    return (mp.e**(x*mp.log(x)) - 1)/x

def dFc(x):
    return x**(x-1)*(mp.log(x) + 1 - 1/x) + 1/x**2

# dominant conjugate-pair fit:  b_n ~ A lam^n + cc(lam)^n  =>  b_{n+2} = S b_{n+1} - P b_n
def fit_pair(k_lo, k_hi):
    """least-squares fit of (S,P) using equations b_{n+2}=S b_{n+1}-P b_n for n in [k_lo,k_hi]"""
    rows = []
    for n in range(k_lo, k_hi):
        rows.append((float(ds[n+1]), -float(ds[n]), float(ds[n+2])))
    # normal equations
    s11 = sum(r[0]*r[0] for r in rows); s12 = sum(r[0]*r[1] for r in rows)
    s22 = sum(r[1]*r[1] for r in rows)
    t1 = sum(r[0]*r[2] for r in rows);  t2 = sum(r[1]*r[2] for r in rows)
    det = s11*s22 - s12*s12
    S = (t1*s22 - t2*s12)/det
    P = (s11*t2 - s12*t1)/det
    return S, P

Sfit, Pfit = fit_pair(45, 58)
disc = Sfit**2 - 4*Pfit
if disc < 0:
    lam = complex(Sfit/2, math.sqrt(-disc)/2)
else:
    lam = complex(Sfit/2 + math.sqrt(disc)/2, 0)
y1 = 1/lam
print(f"  fitted dominant singularity pair: lambda = {lam}")
print(f"  y1 = 1/lambda = {y1}   |y1| = {abs(y1):.8f}   arg = {math.degrees((cmath.phase(y1)))%360:.3f} deg")
OUT['singularity_fit'] = {'lam_real': lam.real, 'lam_imag': lam.imag,
                          'y1_real': y1.real, 'y1_imag': y1.imag, 'R': abs(y1)}
R_fit = abs(y1)
ang_fit = cmath.phase(y1)

# path continuation along the fitted direction, then locate the critical point
def newton_F(z, x0, tol=mp.mpf(10)**-24, steps=60):
    x = x0
    for it in range(steps):
        try:
            dx = (Fc(x) - z)/dFc(x)
        except (ZeroDivisionError, ValueError):
            return x, False
        if not (abs(dx) < mp.mpf(10)**40):
            return x, False
        x = x - dx
        if abs(dx) < tol:
            return x, True
    return x, abs(Fc(x) - z) < tol

ct, st = mp.cos(mp.mpf(ang_fit)), mp.sin(mp.mpf(ang_fit))
xprev = mp.mpc(1)
tt = mp.mpf('0.005')
tmax = R_fit*mp.mpf('0.985')
while tt < tmax:
    z = tt*(ct + 1j*st)
    xr, ok = newton_F(z, xprev)
    if not ok:
        break
    xprev = xr
    tt *= mp.mpf('1.005')
print(f"  path endpoint near singularity: x_end = {xprev}")
# Newton on G from perturbations around x_end and mirrored
Gfun = lambda x: x**(x+1)*(mp.log(x) + 1 - 1/x) + 1
crit = []
seeds = [xprev, xprev*1.05, xprev*0.95, xprev + 0.1, xprev - 0.1,
         xprev + 0.1j, xprev - 0.1j, mp.conj(xprev), mp.conj(xprev)*1.05,
         mp.conj(xprev) + 0.1, mp.conj(xprev) - 0.1j, xprev + 0.05 - 0.07j,
         xprev*1.1, xprev*0.9, xprev*(1+0.1j)]
for seed in seeds:
    for s2 in (seed, seed*1.04):
        try:
            with mp.workdps(50):
                r = mp.findroot(Gfun, (seed, s2), tol=mp.mpf(10)**-30, maxsteps=50)
            if abs(Gfun(r)) < mp.mpf(10)**-20 and abs(r) < 50:
                if all(abs(r - rr) > mp.mpf(10)**-8 for rr, _ in crit):
                    crit.append((r, Fc(r)))
        except Exception:
            pass
print(f"  {len(crit)} distinct critical points found near the path")
crit_s = sorted(crit, key=lambda p: abs(p[1]))
OUT['critical_points'] = [(mp.nstr(z[0], 14), mp.nstr(z[1], 14), float(abs(z[1]))) for z in crit_s[:8]]
for rc, yc in crit_s[:8]:
    print(f"    x_c = {mp.nstr(rc,14)}   y_c = {mp.nstr(yc,14)}   |y_c| = {float(abs(yc)):.10f}")
if crit_s:
    R_crit = float(abs(crit_s[0][1]))
    print(f"  ==> R (nearest critical value) = {R_crit:.10f}   [coefficient fit: {R_fit:.8f}]")
    OUT['R_final'] = R_crit
else:
    OUT['R_final'] = R_fit

# broad sweep for critical points to confirm this pair is the nearest
crit2 = list(crit)
for re_s in [0.25, 0.45, 0.65, 0.85, 1.05, 1.35, 1.75, 2.3]:
    for im_s in [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.9, 2.4]:
        for sgn in (1, -1):
            seed = complex(re_s, sgn*im_s)
            try:
                with mp.workdps(40):
                    r = mp.findroot(Gfun, (seed, seed*1.03), tol=mp.mpf(10)**-25, maxsteps=30)
                if abs(Gfun(r)) < mp.mpf(10)**-16 and abs(r) < 50:
                    if all(abs(r - rr) > mp.mpf(10)**-8 for rr, _ in crit2):
                        crit2.append((r, Fc(r)))
            except Exception:
                pass
crit2_s = sorted(crit2, key=lambda p: abs(p[1]))
print(f"  broad sweep: {len(crit2_s)} distinct critical points total")
OUT['critical_points_broad'] = [(mp.nstr(z[0], 12), mp.nstr(z[1], 12), float(abs(z[1])))
                                 for z in crit2_s[:10]]
for rc, yc in crit2_s[:10]:
    print(f"    x_c = {mp.nstr(rc,12)}   y_c = {mp.nstr(yc,12)}   |y_c| = {float(abs(yc)):.8f}")

# empirical radius check: partial sums inside vs outside R
if crit_s:
    Rv = abs(crit_s[0][1])
    for yv, lab in [(Rv*0.9, 'inside'), (Rv*1.06, 'outside')]:
        partials = []
        acc = mp.mpf(1)
        s = mp.mpf(1)
        for k in range(1, 25):
            s += mp.mpf(ds[k].numerator)/mp.mpf(ds[k].denominator)*mp.mpf(yv)**k
            partials.append(abs(s - T(mp.mpf(float(yv)))))
        print(f"  partials at y={mp.nstr(yv,6)} ({lab}): err[8]={mp.nstr(partials[7],3)}, err[16]={mp.nstr(partials[15],3)}, err[24]={mp.nstr(partials[23],3)}")

# ---------------------------------------------------------------- monotonicity / sandwich numerics
print("\n== Theorem numerics ==")
mono_ok = True
for i in range(1, 240):
    x = mp.mpf(10)**(mp.mpf(-6) + (mp.mpf(2) - mp.mpf(-6)) * i / 240)
    fp = x**(x-1)*(mp.log(x) + 1 - 1/x) + 1/x**2
    if fp <= 0:
        mono_ok = False
print("f'(x) > 0 on grid (1e-6,100):", mono_ok)
OUT['monotone_grid_ok'] = mono_ok

sand1 = all(T(y) < mp.e**y for y in [mp.mpf(k)/8 for k in range(-80, 80) if k != 0])
sand2 = all(1 + y < T(y) for y in [mp.mpf(k)/40 for k in range(-39, 0)])
sand3 = all(T(y) < 1 + y for y in [mp.mpf(k)/40 for k in range(1, 400)])
print("T(y)<e^y (y != 0):", sand1, "| 1+y<T(y) on (-1,0):", sand2, "| T(y)<1+y on (0,10):", sand3)
OUT['sandwich_ok'] = bool(sand1 and sand2 and sand3)

# T' <= 1 check: h(x) = x^x[x(lnx+1)-1]+1-x^2 >= 0
h_ok = True
for i in range(1, 200):
    x = mp.mpf(10)**(mp.mpf(-4) + (mp.mpf(2) - mp.mpf(-4)) * i / 200)
    h = x**x * (x*mp.log(x) + x - 1) + 1 - x**2
    if h < -mp.mpf(10)**-30:
        h_ok = False
print("f'(x) >= 1 grid check (h >= 0):", h_ok)
OUT['lipschitz_grid_ok'] = h_ok

# ---------------------------------------------------------------- golden comparison table
print("\n== Golden/metallic comparison ==")
rows = []
for b in ['-1', '-0.5', '0.5', '1', '1.4', '1.49', '1.5', '1.51', '2', '3', '10']:
    bb = mp.mpf(b)
    Tb, ph = T(bb), phi_metal(bb)
    rows.append((float(bb), float(Tb), float(ph), float(Tb - ph)))
    print(f"  b={b:>5}: T(b)={mp.nstr(Tb,10):>12}  phi_b={mp.nstr(ph,10):>10}  diff={mp.nstr(Tb-ph,4):>10}")
OUT['golden_table'] = rows
print("exact: T(0)=1=phi_0 ->", T(0) == 1, " T(3/2)=2=phi_{3/2} ->", T(mp.mpf('1.5')) == 2,
      "| phi_{3/2} =", mp.nstr(phi_metal(mp.mpf('1.5')), 15))

# ---------------------------------------------------------------- asymptotics
print("\n== Asymptotics ==")
asym = {}
for Y in [-5, -10, -20, -40]:
    Y = mp.mpf(Y)
    x = T(Y)
    asym[f'y={Y}'] = [mp.nstr(x / mp.e**Y, 12), mp.nstr(mp.log(x/mp.e**Y), 8)]
    print(f"  y={mp.nstr(Y,4)}: T/e^y = {asym[f'y={Y}'][0]},  ln(T/e^y) = {asym[f'y={Y}'][1]}   (pred -y^2 e^y/2 = {mp.nstr(-Y**2*mp.e**Y/2,6)})")
OUT['asym_neg'] = asym

def ssqrt(yv):
    return mp.e**mp.lambertw(mp.log(yv))
asym2 = {}
for Y in [10**3, 10**6, 10**9, 10**12, 10**18]:
    Y = mp.mpf(Y)
    x = T(Y)
    s = ssqrt(Y)
    asym2[str(Y)] = [mp.nstr(x, 10), mp.nstr(x - s, 10), mp.nstr(x*mp.log(mp.log(Y))/mp.log(Y), 8)]
    print(f"  y=10^{math.log10(float(Y))}: T={mp.nstr(x,8)},  T-ssqrt={mp.nstr(x-s,8)},  T*lnln/ln={mp.nstr(x*mp.log(mp.log(Y))/mp.log(Y),6)}")
OUT['asym_pos'] = asym2
print("  identity x^(x-1) = y + 1/x check:",
      mp.nstr(T(mp.mpf(100))**(T(mp.mpf(100))-1) - mp.mpf(100) - 1/T(mp.mpf(100)), 5))
print("  T(10^100) =", mp.nstr(T(mp.mpf(10)**100), 12))
OUT['T_googol'] = mp.nstr(T(mp.mpf(10)**100), 12)
print("  T(10^6)   =", mp.nstr(T(mp.mpf(10)**6), 12))

# ---------------------------------------------------------------- exact special values (sympy exact)
print("\n== Exact special values (symbolic verification) ==")
xs = sp.symbols('x', positive=True)
def fs(z):  # sympy f at rational/surd z
    return sp.simplify((z**z - 1)/z)
checks = [
    ("T(0)=1",        sp.Rational(0), sp.Integer(1)),
    ("T(3/2)=2",      sp.Rational(3,2), sp.Integer(2)),
    ("T(26/3)=3",     sp.Rational(26,3), sp.Integer(3)),
    ("T(255/4)=4",    sp.Rational(255,4), sp.Integer(4)),
    ("T(sqrt2-2)=1/2", sp.sqrt(2)-2, sp.Rational(1,2)),
    ("T(2sqrt2-4)=1/4", 2*sp.sqrt(2)-4, sp.Rational(1,4)),
    ("T(3^(2/3)-3)=1/3", 3**(sp.Rational(2,3))-3, sp.Rational(1,3)),
]
for label, b, xval in checks:
    lhs = sp.simplify(fs(xval) - b)
    print(f"  {label}: {sp.simplify(lhs) == 0}")
    OUT[f'exact_{label}'] = bool(sp.simplify(lhs) == 0)

# negative integer roots: x=-n solves x^x = b_n x + 1,  b_n=(1-(-1)^n n^-n)/n
print("\n== Negative-integer exact roots of the super-quadratic ==")
bn_list = []
for n in range(1, 8):
    bn = sp.Rational(1 - ((-1)**n) * sp.Integer(n)**(-n)) / n
    xneg = sp.Integer(-n)
    ok = sp.simplify(xneg**xneg - bn*xneg - 1) == 0
    bn_list.append((n, str(bn), float(bn), bool(ok)))
    print(f"  n={n}: b_n = {bn} = {float(bn):.12f}   (-n)^(-n)=b_n(-n)+1 : {ok}")
OUT['bn_family'] = bn_list

# ---------------------------------------------------------------- dynamics
print("\n== Dynamics ==")
iters = []
xk = mp.e**mp.lambertw(mp.log(2))
for k in range(14):
    xk = mp.e**mp.lambertw(mp.log(xk + 1))
    iters.append(mp.nstr(xk, 18))
    print(f"  W-iter {k}: x = {mp.nstr(xk,18)}   err = {mp.nstr(xk - alpha, 5)}")
OUT['W_iteration'] = iters
gprime_alpha = (alpha+1)*mp.log(alpha)     # multiplier of x->alpha^x-1 at alpha
OUT['tower_multiplier_alpha'] = mp.nstr(gprime_alpha, 12)
print("  multiplier of g(x)=alpha^x-1 at alpha:", OUT['tower_multiplier_alpha'])
# tower iteration from 0.5 -> 0
tseq = []
xk = mp.mpf('0.5')
for k in range(30):
    xk = alpha**xk - 1
    if k < 12 or k in (15, 20, 29):
        tseq.append((k, mp.nstr(xk, 10)))
print("  tower iter x->alpha^x-1 from 0.5:", tseq[:6], "... ->", tseq[-1])
OUT['tower_to_zero'] = tseq
# orbit of T from 1 -> delta
orb = []
yk = mp.mpf(1)
for k in range(9):
    yk = T(yk)
    orb.append(mp.nstr(yk, 14))
print("  orbit T^n(1):", orb)
OUT['orbit_T'] = orb
fdelta = delta**(delta+1)*(mp.log(delta)+1) - delta**delta + 1
OUT['Tprime_delta'] = mp.nstr(delta**2/fdelta, 12)
print("  T'(delta) =", OUT['Tprime_delta'], " (attracting: <1)")
# Newton for alpha
nw = []
xk = mp.mpf(2)
for k in range(7):
    xk = xk - (xk**xk - xk - 1)/(xk**xk*(mp.log(xk)+1) - 1)
    nw.append(mp.nstr(xk, 22))
print("  Newton from 2:", nw)
OUT['newton_alpha'] = nw

# ---------------------------------------------------------------- integral
print("\n== Integrals ==")
I_direct = mp.quad(lambda yy: T(yy), [0, 1])
I_redux = alpha + mp.log(alpha) - mp.quad(lambda x: x**(x-1), [1, alpha])
print("  int_0^1 T  =", mp.nstr(I_direct, 15))
print("  redux alpha+ln alpha - int x^(x-1) =", mp.nstr(I_redux, 15))
OUT['integral_0_1'] = mp.nstr(I_direct, 15)
OUT['integral_check'] = mp.nstr(I_direct - I_redux, 5)
I_neg = mp.quad(lambda yy: 1 - T(yy), [-1, 0])
OUT['integral_neg'] = mp.nstr(I_neg, 15)
print("  int_{-1}^0 (1-T) =", OUT['integral_neg'])

# ---------------------------------------------------------------- numerology
print("\n== Numerology (near misses) ==")
num = {
    'alpha_vs_sqrtpi': mp.nstr(alpha - mp.sqrt(mp.pi), 12),
    'alpha_plus_inv_vs': mp.nstr(alpha + 1/alpha - (mp.sqrt(mp.pi) + mp.lambertw(1)), 12),
    'alpha_beta_prod': mp.nstr(alpha*beta, 12),
    'alpha_beta_sum': mp.nstr(alpha + beta, 12),
    'alpha_minus_phi': mp.nstr(alpha - phi, 12),
    'inv_alpha': mp.nstr(1/alpha, 12),
}
for k, s in num.items():
    print(f"  {k} = {s}")
OUT['numerology'] = num

# ---------------------------------------------------------------- analytic lemma grid checks
print("\n== Analytic convexity/Lipschitz lemma checks (grid) ==")
l_ok = all(mp.e**(x*mp.log(x)) <= 1/(1 - x*mp.log(x)) for x in [mp.mpf(k)/500 for k in range(1, 500)])
r_ok = all(1 - x + x**2*mp.log(x) > 0 for x in [mp.mpf(k)/500 for k in range(1, 500)])
u_ok = all(x*(mp.log(x)+1)**2 < 1 for x in [mp.mpf(k)/500 for k in range(1, 500)])
p_ok = all(x*mp.log(x)*(mp.log(x)+1) - x + 1 >= 0 for x in [1 + mp.mpf(k)/300 for k in range(0, 400)])
Pgrid_ok = all(x**x*((x*mp.log(x)+x-1)**2 + x + 1) < 2 for x in [mp.mpf(k)/600 for k in range(1, 600)])
hgrid_ok = all(x**x*(x*mp.log(x)+x-1) + 1 - x**2 >= 0 for x in [mp.mpf(10)**(mp.mpf(-5) + mp.mpf(2)*k/300) for k in range(1, 300)])
print("  e^t <= 1/(1-t) chain:", l_ok, "| r(x)=1-x+x^2 lnx > 0:", r_ok,
      "| x(lnx+1)^2 < 1:", u_ok)
print("  p(x) >= 0 on [1,2.33]:", p_ok, "| P(x) < 2 grid:", Pgrid_ok,
      "| h(x) >= 0 grid:", hgrid_ok)
OUT['lemma_checks'] = bool(l_ok and r_ok and u_ok and p_ok and Pgrid_ok and hgrid_ok)

# ---------------------------------------------------------------- value table
print("\n== Value table ==")
table = []
for Y in ['-10', '-5', '-2', '-3/2', '-1', '-1/2', '0', '1/2', '1', '3/2', '2', 'e', 'pi', '5', '10', '100', '10^6']:
    Yv = mp.mpf(Y.replace('10^6', '1000000').replace('e', str(mp.e)).replace('pi', str(mp.pi)))
    val = T(Yv)
    row = [Y, mp.nstr(val, 18)]
    table.append(row)
    print(f"  T({Y:>6}) = {mp.nstr(val, 18)}")
OUT['value_table'] = table

with open('tahawus_results.json', 'w') as fh:
    json.dump(OUT, fh, indent=1, default=str)
print("\nSaved -> tahawus_results.json")
