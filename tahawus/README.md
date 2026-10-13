# The Tahawus Function

**An Extension of the Lambert W Function**

*Justo Tapiador García — Universidad de Alicante (UA), Spain*
Preprint, September 27, 2026 · 29 pp. + cover · [PDF](tahawus.pdf)
MSC 2020: 26A18 (primary); 33E30, 11J81, 37C25, 65H05.

An expanded and rigorous version of the author's 2018 blog post
*"Hola, os presento a Tahawus: una extensión de la función W de Lambert"*
([Tardigrados](https://tardigrados.wordpress.com/2018/04/03/hola-os-presento-a-tahawus-una-extension-de-la-funcion-w-de-lambert)).
The function is named after Tahawus, New York — the Adirondack ghost town
where Richard P. Stanley grew up.

## Abstract

The *Tahawus function* T is the inverse of the map f(x) = (x^x − 1)/x. It
was introduced in the 2018 post as a natural companion of the Lambert W
function: T solves the "super-quadratic" x^x = bx + 1 the same way W solves
x·e^x = a. This paper gives the function a systematic theory.

## Main results

- **Analytic theory.** f is a strictly increasing real-analytic bijection
  (0, ∞) → R, so T: R → (0, ∞) is a real-analytic order isomorphism.
- **Sandwich theorem.** 1 + y < T(y) < e^y on (−1, 0), while
  T(y) < min(1+y, e^y) on (0, ∞); with T(y) ~ e^y as y → −∞,
  T(y) ~ ln y / ln ln y as y → +∞, and the refined law
  T(y) − e^{W(ln y)} → 1.
- **Arithmetic.** For every rational b, T(b) is either a positive integer
  (exactly when b = (n^n − 1)/n) or a transcendental number. This recovers
  Stanley's 1989 observation that the *Tahawus constant*
  α = T(1) = 1.7767750400…, the unique real root of x^x = x + 1, is
  transcendental.
- **The constellation.** Five connected constants, all in the OEIS, all
  transcendental: α = T(1) (A124930), β = T(−1) with β^β + β = 1
  (A226568), c = 1/β = 3.2931… (A169862, c^(c−1) = (c−1)^c),
  v = c − 1 = 2.2931… (A085846; v and v+1 are the unique pair of
  consecutive positive reals with equal x^(1/x)), and υ₂ = 5.5079…
  (A144184, root of x^(1+1/x) = x + 2).
- **Local theory.** The exact Maclaurin expansion of T (whose quadratic
  term vanishes), a convexity dichotomy, and the numerical determination of
  the radius of convergence R ≈ 0.65261.
- **Golden comparison.** T(b) versus the metallic-ratio family
  φ_b = (b + √(b²+4))/2, with the unexpected exact coincidences
  T(0) = φ_0 = 1 and T(3/2) = φ_{3/2} = 2; plus the negative-integer exact
  roots x = −n of the super-quadratic family.
- **Dynamics.** A unique attracting fixed point δ = 2.2288… (root of
  x^x = x² + 1, transcendental), a convergent reading of the blog post's
  nested composition, and a precise correction of the original post's
  infinite-tower identity.

All numerical claims are verified to 70 significant digits.

## Contents

| File | Description |
|------|-------------|
| `tahawus.pdf` | The complete paper: designed cover + body (30 pp.), with PDF bookmarks. |
| `main.pdf` | The paper body alone (29 pp.): table of contents, title block, sections, bibliography. |
| `cover.pdf` | The designed cover page. |
| `main.tex` + `sections/` | LaTeX source: 7 section files covering the 11 sections and 2 appendices. Compile from this folder with `tectonic main.tex`. |
| `figures/` | The six vector figures (mirror, sandwich, growth, slit, orbit, constants). |
| `scripts/tahawus_research.py` | Produces all numerical facts of the paper at 80-digit precision; writes `tahawus_results.json`. |
| `scripts/tahawus_extra.py` | Additional numerical verifications at 50 digits; writes `tahawus_extra.json`. |
| `scripts/tahawus_figures.py` | Regenerates the six figures (run from `scripts/`; writes to `../figures`). |
| `scripts/*.json` | The recorded numerical results (57 + 17 entries). |

## Reproduction

```bash
pip install mpmath sympy numpy matplotlib
cd scripts
python3 tahawus_research.py   # -> tahawus_results.json
python3 tahawus_extra.py      # -> tahawus_extra.json
python3 tahawus_figures.py    # -> ../figures/*.pdf
```

The paper itself recompiles with `tectonic main.tex` (or `pdflatex`, two
passes) from this folder; the figures must remain in `figures/`.

## Citation

BibTeX:

```bibtex
@unpublished{tapiador2026tahawus,
  author = {Tapiador Garc{\'i}a, Justo},
  title  = {The Tahawus Function: An Extension of the Lambert W Function},
  year   = {2026},
  note   = {Preprint, available at
            https://github.com/Justo-Tapiador/maths/blob/main/tahawus/tahawus.pdf}
}
```

## License

Code and data: MIT (see the repository LICENSE). Paper:
© Justo Tapiador García.
