# The Prime Staircase

**Superexponential Growth along the Diagonal of the Hierarchy of
Higher-Order Primes**

*Justo Tapiador — Universidad de Alicante (UA), Spain*
Preprint, October 2026 · 12 pages · [PDF](prime_staircase.pdf)
MSC 2020: 11A41 (primary); 11N05, 11Y05.

## Abstract

Let p_n denote the n-th prime number and let f(n) = p_n be the
prime-indexing (primeth) map. Iterating the selection of prime-indexed terms
generates the classical hierarchy of higher-order primes
p_{n,k} = f^{k+1}(n), whose first levels are the primes themselves, the
superprimes of Dressler and Parker (A006450), and their higher-order
analogues (A038580, A049090, ...). This paper studies the diagonal of that
hierarchy, d_k = p_{k,k} = f^{k+1}(k), obtained by taking from each level k
the single term whose index coincides with the level — named here the
**prime staircase**:

```
3, 11, 127, 1787, 52711, 1128889, 37139213, 997525853,
36294260117, 1765037224331, 75063692618249, 4176603711876241,
214616508361122337, 10643741922535940459
```

The sequence was not catalogued in the OEIS when the paper was written
(searches of October 2026); a submission is in preparation.

## Main results

- **Master formula.** li(log d_k) = k + 1 + li(log k) + O(log log k), with
  numerical residuals below 0.2 for k ≥ 4. This yields the second-order
  expansion log d_k = k(log k + log log k − 1) + (log k)(1 + li(log k)) +
  O(k(log log k)²/log k).
- **Growth.** In particular log d_k ~ k log k ~ p_k: the logarithm of the
  k-th staircase prime is asymptotic to the k-th prime itself. For k ≥ 11,
  d_k has exactly 2k − 8 decimal digits.
- **Structure.** The order-of-primeness ladder F(d_k) = k + 1 + F(k), where
  F = A049076 is the order of primeness; and the orbit-crossing law: d_k
  lies on the primeth recurrence A007097 if and only if k does (first
  crossings at k = 31, 127, 709, ...). The full iteration matrix is
  A354967, of which this sequence is the superdiagonal; its index twin
  m_k = f^k(k) is A058009, studied by Miska, Tóth and Żmija
  (arXiv:1908.10421).
- **Computation.** d_1, ..., d_10 are computed exactly in the paper (the
  largest, d_10 = 1765037224331, in under four minutes on two cores), and
  d_11, ..., d_14 are bracketed rigorously there. The exact values of
  d_11, ..., d_14 listed above were subsequently confirmed by independent
  prime-counting computations (see Verification), which situates the
  exact-computability horizon of the sequence around k = 14 with current
  π(x) technology.

## Contents

| File | Description |
|------|-------------|
| `prime_staircase.pdf` | The paper (12 pp., vector figures). |
| `prime_staircase.tex` | LaTeX source. Compile with `tectonic prime_staircase.tex` (or `pdflatex`, 2 passes); the figure PDFs must sit next to it. |
| `fig_growth.pdf` | Figure 1: growth of the staircase. |
| `fig_logd_pk.pdf` | Figure 2: log d_k versus the k-th prime. |
| `staircase_verify.py` | Self-contained verification script (requires only `sympy`): recomputes the hierarchy and its OEIS levels, m_k against A058009, the staircase d_1..d_8, the structural propositions and the master-formula residuals. Runs in a few seconds. |
| `staircase_data.json` | Numerical data: hierarchy rows, chains, m_k, d_k, growth tables, and the rigorous brackets for d_11..d_14. |

## Verification

```bash
pip install sympy
python3 staircase_verify.py
```

This recomputes d_1..d_8 from scratch (d_9 and d_10, obtained by segmented
sieving, are reproduced as constants) and checks the structural and
asymptotic results numerically. The heavy computations behind d_9..d_14 use
the `primesieve` and `primecount` libraries of Kim Walisch, documented in
Section 5 of the paper; the exact values of d_11..d_14 were confirmed by
four independent checks: prime-counting cross-validation against the
published values of pi(10^13) and pi(10^14), full orbit rebuilds landing on
A058009(11..14), and consistency with the rigorous brackets recorded in
`staircase_data.json`.

## Citation

BibTeX:

```bibtex
@unpublished{tapiador2026staircase,
  author = {Tapiador, Justo},
  title  = {The Prime Staircase: Superexponential Growth along the Diagonal
            of the Hierarchy of Higher-Order Primes},
  year   = {2026},
  note   = {Preprint, available at
            https://github.com/Justo-Tapiador/maths/blob/main/prime-staircase/prime_staircase.pdf}
}
```

## License

Code and data: MIT (see the repository LICENSE). Paper: © Justo Tapiador.
