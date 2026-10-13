#!/usr/bin/env python3
"""Figures for the Tahawus function paper. Vector PDF output, English labels."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import mpmath as mp

for fp in ['/usr/share/fonts/truetype/english/Tinos-Regular.ttf',
           '/usr/share/fonts/truetype/english/Tinos-Bold.ttf',
           '/usr/share/fonts/truetype/english/Tinos-Italic.ttf']:
    try:
        fm.fontManager.addfont(fp)
    except Exception:
        pass
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Tinos', 'Liberation Serif', 'DejaVu Serif'],
    'mathtext.fontset': 'stix',
    'font.size': 10.5,
    'axes.labelsize': 11,
    'axes.linewidth': 0.8,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': False, 'ytick.right': False,
    'legend.frameon': False,
    'figure.dpi': 150,
})
mp.mp.dps = 30
OUT = '../figures'

C_F, C_T, C_E, C_L, C_G, C_O = '#1f5fa8', '#c0392b', '#7d3c98', '#5d6d7e', '#1e8449', '#ca6f1e'

def f(x):
    return (x**x - 1)/x

def T(yv):
    y = mp.mpf(yv)
    if y > 0:
        lo, hi = mp.mpf(1), mp.mpf(2)
        while f(hi) < y:
            hi *= 2
    else:
        lo, hi = mp.mpf('1e-25'), mp.mpf(1)
    for _ in range(120):
        mid = (lo + hi)/2
        if f(mid) < y:
            lo = mid
        else:
            hi = mid
    try:
        return mp.findroot(lambda x: x**x - y*x - 1, (lo, hi), tol=mp.mpf(10)**-25, maxsteps=40)
    except Exception:
        return (lo + hi)/2

alpha = T(1); beta = T(-1); delta = mp.findroot(lambda x: x**x - x**2 - 1, (mp.mpf(2.2), mp.mpf(2.3)))

def style(ax):
    ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# ---------------------------------------------------------------- Figure 1: mirror
fig, ax = plt.subplots(figsize=(6.4, 5.4), constrained_layout=True)
xs = np.linspace(0.045, 2.75, 500)
fs = [float(f(mp.mpf(x))) for x in xs]
ax.plot(xs, fs, color=C_F, lw=2.2, label=r'$f(x)=\dfrac{x^{x}-1}{x}$')
ys = np.linspace(-3.2, 3.9, 700)
Ts = [float(T(mp.mpf(y))) for y in ys]
ax.plot(ys, Ts, color=C_T, lw=2.2, label=r'$x=\mathcal{T}(y)$')
ax.plot([-1, 4], [-1, 4], color=C_L, lw=1, ls='--', alpha=0.8)
ax.plot([0, 2.2], [-1, 1.2], color=C_L, lw=0.9, ls=':', alpha=0.9)   # y = x - 1
ax.plot([-1, 2.2], [0, 3.2], color=C_L, lw=0.9, ls=':', alpha=0.9)   # y = x + 1
ax.axhline(0, color='gray', lw=0.6); ax.axvline(0, color='gray', lw=0.6)
for pt, lab, dx, dy in [((1, 0), r'$(1,0)$', 8, -14), ((0, 1), r'$(0,1)$', 6, 4),
                        ((alpha, 1), r'$(\alpha,1)$', 6, 4), ((1, alpha), r'$(1,\alpha)$', -20, 2),
                        ((-1, beta), r'$(-1,\beta)$', -4, 8), ((delta, delta), r'$(\delta,\delta)$', 8, -4)]:
    ax.plot(*pt, 'o', color='black', ms=4.5)
    ax.annotate(lab, pt, textcoords='offset points', xytext=(dx, dy), fontsize=10)
ax.set_xlim(-1.2, 4); ax.set_ylim(-3.2, 4)
ax.set_xlabel('$x$'); ax.set_ylabel('$y$')
ax.legend(loc='lower right', fontsize=11)
style(ax)
fig.savefig(f'{OUT}/fig1_mirror.pdf'); plt.close(fig)
print('fig1 done')

# ---------------------------------------------------------------- Figure 2: sandwich
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 4.1), constrained_layout=True)
ys = np.linspace(-4, 2, 400)
Ts = np.array([float(T(mp.mpf(y))) for y in ys])
a1.plot(ys, np.exp(ys), color=C_E, lw=2, label=r'$e^{y}$')
a1.plot(ys, Ts, color=C_T, lw=2.4, label=r'$\mathcal{T}(y)$')
a1.plot(ys, 1 + ys, color=C_G, lw=2, label=r'$1+y$')
a1.fill_between(ys, 1 + ys, np.exp(ys), color=C_G, alpha=0.07)
a1.axhline(1, color='gray', lw=0.5, ls=':'); a1.axvline(0, color='gray', lw=0.5, ls=':')
a1.plot([0], [1], 'o', color='black', ms=4.5)
a1.annotate(r'$\mathcal{T}(0)=1$', (0, 1), textcoords='offset points', xytext=(10, -14), fontsize=10)
a1.set_xlim(-4, 2); a1.set_ylim(-0.05, 3.6)
a1.set_xlabel('$y$'); a1.set_title(r'(a) $\mathcal{T}$ squeezed between $1+y$ and $e^y$', fontsize=10.5)
a1.legend(loc='upper left', fontsize=10)
style(a1)
ys2 = np.linspace(-30, -1.0, 500)
rat = np.array([float(T(mp.mpf(y))/mp.e**mp.mpf(y)) for y in ys2])
a2.plot(ys2, rat, color=C_T, lw=2.4, label=r'$\mathcal{T}(y)/e^{y}$')
a2.axhline(1, color=C_L, lw=1, ls='--')
asym = 1 - (ys2**2/2)*np.exp(ys2)
a2.plot(ys2, asym, color=C_E, lw=1.4, ls=':', label=r'$1-\frac{y^{2}}{2}e^{y}$')
a2.set_xlabel('$y$'); a2.set_xlim(-30, 0)
a2.set_title(r'(b) approach to $e^y$ as $y\to-\infty$', fontsize=10.5)
a2.legend(loc='upper left', fontsize=10)
style(a2)
fig.savefig(f'{OUT}/fig2_sandwich.pdf'); plt.close(fig)
print('fig2 done')

# ---------------------------------------------------------------- Figure 3: growth at +infty
fig, ax = plt.subplots(figsize=(6.4, 4.6), constrained_layout=True)
yy = np.logspace(1, 9, 160)
Tv = np.array([float(T(mp.mpf(y))) for y in yy])
ss = np.array([float(mp.e**mp.lambertw(mp.log(mp.mpf(y)))) for y in yy])
ax.loglog(yy, Tv, color=C_T, lw=2.4, label=r'$\mathcal{T}(y)$')
ax.loglog(yy, ss, color=C_E, lw=1.8, ls='--', label=r'$e^{W(\ln y)}$ (super square root)')
ax.loglog(yy, np.log(yy), color=C_G, lw=1.8, label=r'$\ln y$')
ax.loglog(yy, np.log(yy)/np.log(np.log(yy)), color=C_L, lw=1.6, ls=':', label=r'$\ln y/\ln\ln y$')
ax.set_xlabel('$y$'); ax.set_ylabel(r'$\mathcal{T}(y)$')
ax.set_xlim(10, 1e9); ax.set_ylim(1.8, 30)
ax.legend(loc='upper left', fontsize=10)
ax.set_title(r'$\mathcal{T}$ grows like $\ln y/\ln\ln y$; $\mathcal{T}(y)\sim e^{W(\ln y)}+1$', fontsize=10.5)
style(ax)
fig.savefig(f'{OUT}/fig3_growth.pdf'); plt.close(fig)
print('fig3 done')

# ---------------------------------------------------------------- Figure 4: image of the slit
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 4.3), constrained_layout=True)
r = np.linspace(0.28, 40, 3000)
gam = np.array([complex((1 - (mp.mpf(rr)**(-mp.mpf(rr)))*mp.e**(-1j*mp.pi*mp.mpf(rr)))/mp.mpf(rr))
                for rr in r])
a1.plot(gam.real, gam.imag, color=C_E, lw=1.8, label=r'$\gamma^{+}: F(-r)$, $r>0$')
a1.plot(gam.real, -gam.imag, color=C_E, lw=1.2, ls='--', alpha=0.65, label=r'$\overline{\gamma^{+}}$')
th = np.linspace(0, 2*np.pi, 200)
a1.plot(0.6526*np.cos(th), 0.6526*np.sin(th), color=C_L, lw=1, ls=':', alpha=0.9)
a1.plot(-0.104895, 0.644125, 's', color=C_T, ms=6)
a1.plot(-0.104895, -0.644125, 's', color=C_T, ms=6)
a1.annotate(r'$y_c,\ \bar y_c$', (-0.104895, 0.644125), textcoords='offset points',
            xytext=(8, 6), fontsize=10, color=C_T)
a1.annotate(r'$|y|=R$', (0.6526*np.cos(0.8), 0.6526*np.sin(0.8)), textcoords='offset points',
            xytext=(2, 8), fontsize=9.5, color=C_L)
for n in range(1, 9):
    bn = (1 - ((-1)**n)*n**(-n))/n
    a1.plot(bn, 0, 'o', color='black', ms=3.6)
a1.annotate(r'$b_n$', ((1-0.5**-1)/2*0.5+0.42, 0.05), fontsize=10)
a1.plot(0, 0, 'x', color='black', ms=6)
a1.set_xlabel(r'$\Re y$'); a1.set_ylabel(r'$\Im y$')
a1.set_xlim(-0.6, 2.3); a1.set_ylim(-1.7, 1.7)
a1.axhline(0, color='gray', lw=0.5); a1.axvline(0, color='gray', lw=0.5)
a1.set_title(r'(a) image of the slit winds into $0$', fontsize=10.5)
a1.legend(loc='upper right', fontsize=9.5)
style(a1)
a2.plot(gam.real, gam.imag, color=C_E, lw=1.6)
a2.plot(gam.real, -gam.imag, color=C_E, lw=1.1, ls='--', alpha=0.65)
a2.plot(0.6526*np.cos(th), 0.6526*np.sin(th), color=C_L, lw=1.1, ls=':')
a2.plot(-0.104895, 0.644125, 's', color=C_T, ms=6)
a2.plot(-0.104895, -0.644125, 's', color=C_T, ms=6)
for n in range(2, 12):
    bn = (1 - ((-1)**n)*n**(-n))/n
    a2.plot(bn, 0, 'o', color='black', ms=3.4)
a2.plot(0, 0, 'x', color='black', ms=7)
a2.annotate(r'$b_n\to 0$', (0.12, 0.1), fontsize=10)
a2.set_xlim(-0.06, 0.55); a2.set_ylim(-0.75, 0.75)
a2.set_xlabel(r'$\Re y$')
a2.set_title(r'(b) zoom: $b_n$ accumulate at the origin; $|y|=R$', fontsize=10.5)
style(a2)
fig.savefig(f'{OUT}/fig4_slit.pdf'); plt.close(fig)
print('fig4 done')

# ---------------------------------------------------------------- Figure 5: dynamics
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 4.2), constrained_layout=True)
yy = np.linspace(-1, 3.3, 300)
Tv = np.array([float(T(mp.mpf(y))) for y in yy])
a1.plot(yy, Tv, color=C_T, lw=2.2)
a1.plot(yy, yy, color=C_L, lw=1, ls='--')
yk = 1.0
xs_c, ys_c = [yk], [0]
for k in range(6):
    tn = float(T(mp.mpf(yk)))
    xs_c += [yk, tn]; ys_c += [tn, tn]
    yk = tn
a1.plot(xs_c, ys_c, color=C_O, lw=1.3, alpha=0.9)
for k, (xx, yyv) in enumerate(zip(xs_c[::2], ys_c[::2])):
    a1.plot(xx, yyv, 'o', color=C_O, ms=3.8)
a1.plot(delta, delta, 'o', color='black', ms=5)
a1.annotate(r'$\delta\approx 2.2288$', (delta, delta), textcoords='offset points',
            xytext=(-88, 8), fontsize=10.5)
a1.annotate(r'$y_0=1$', (1, 0), textcoords='offset points', xytext=(-4, -14), fontsize=10)
a1.set_xlabel('$y$'); a1.set_ylabel(r'$\mathcal{T}(y)$')
a1.set_xlim(-1, 3.3); a1.set_ylim(0, 3.3)
a1.set_title(r'(a) every orbit of $\mathcal{T}$ converges to $\delta$', fontsize=10.5)
style(a1)
yy2 = np.linspace(-1, 3.6, 400)
diff = np.array([float(T(mp.mpf(y)) - y) for y in yy2])
a2.plot(yy2, diff, color=C_G, lw=2.2, label=r'$\mathcal{T}(y)-y$')
a2.axhline(0, color='gray', lw=0.6)
a2.axvline(delta, color=C_L, lw=0.9, ls=':')
a2.plot(delta, 0, 'o', color='black', ms=5)
a2.annotate(r'$\delta$: root of $x^{x}=x^{2}+1$', (delta, 0), textcoords='offset points',
            xytext=(-160, 10), fontsize=10)
a2.set_xlabel('$y$'); a2.set_ylabel(r'$\mathcal{T}(y)-y$')
a2.set_title(r'(b) unique fixed point', fontsize=10.5)
a2.legend(loc='lower left', fontsize=10)
style(a2)
fig.savefig(f'{OUT}/fig5_orbit.pdf'); plt.close(fig)
print('fig5 done')

# ---------------------------------------------------------------- Figure 6: three constants
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 4.2), constrained_layout=True)
xx = np.linspace(0.3, 6.0, 500)
hv = np.array([float(x**(1.0/x)) for x in xx])
a1.plot(xx, hv, color=C_E, lw=2.2, label=r'$h(x)=x^{1/x}$')
vv = 1/beta - 1
a1.plot([vv, vv+1], [vv**(1/vv), (vv+1)**(1/(vv+1))], 'o', color='black', ms=5)
a1.annotate(r'$v$', (vv, vv**(1/vv)), textcoords='offset points', xytext=(-16, 2), fontsize=11)
a1.annotate(r'$v+1$', (vv+1, (vv+1)**(1/(vv+1))), textcoords='offset points', xytext=(4, 8), fontsize=11)
a1.plot([0.3, 6], [vv**(1/vv), vv**(1/vv)], color=C_L, lw=0.9, ls=':')
a1.set_xlabel('$x$'); a1.set_ylabel(r'$x^{1/x}$')
a1.set_title(r'(a) $v$ and $v+1$ share the value of $x^{1/x}$', fontsize=10.5)
a1.legend(loc='upper right', fontsize=10)
style(a1)
xx2 = np.linspace(1.0, 3.3, 500)
a2.plot(xx2, np.array([x**x for x in xx2]), color=C_T, lw=2.2, label=r'$x^{x}$')
a2.plot(xx2, np.array([x**(1+1.0/x) for x in xx2]), color=C_E, lw=2.2, label=r'$x^{1+1/x}$')
a2.plot(xx2, xx2 + 1, color=C_G, lw=1.8, label=r'$x+1$')
phi = (1+np.sqrt(5))/2
for pt, lab, off in [(phi, r'$\varphi$', (-18, 22)), (alpha, r'$\alpha$', (-30, -18)), (vv, r'$v$', (6, -16))]:
    a2.plot(pt, pt+1, 'o', color='black', ms=5)
    a2.annotate(lab, (pt, pt+1), textcoords='offset points', xytext=off, fontsize=11)
a2.set_xlabel('$x$'); a2.set_ylabel('$y$')
a2.set_title(r'(b) three constants, one line: $\varphi<\alpha<v$', fontsize=10.5)
a2.legend(loc='upper left', fontsize=10)
style(a2)
fig.savefig(f'{OUT}/fig6_constants.pdf'); plt.close(fig)
print('fig6 done')
print('ALL FIGURES DONE')
