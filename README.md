# A subquadratic upper bound for the mixing time of the Thorp shuffle

This repository contains the manuscript and reproducibility material for

> **A subquadratic upper bound for the mixing time of the Thorp shuffle**  
> Yunjiang Jiang

For a deck of size \(n=2^d\), the paper proves the full-deck worst-start total-variation upper bound

\[
t_{\mathrm{mix}}(\varepsilon)
\le
5.1536951454\ldots\, d^{5/3}
+O_\varepsilon\!\left(d^{4/3}\log(d+2)\right).
\]

Time is counted in individual Thorp shuffles (equivalently, one hypercube sweep consists of \(d\) shuffles).

## What is in this repository

- `paper/main.tex` — self-contained LaTeX source; it is compiled automatically by GitHub Actions, the workflow publishes the PDF as a downloadable build artifact.
- `checks/verify_audit.py` — independent finite-case and numerical audit checks used during proof development.
- `checks/audit_results.json` — machine-readable audit results.
- `checks/run.log` — recorded verifier status summary.
- `notes/SUBMISSION_NOTES.md` — arXiv/submission notes.
- `notes/REVISION_NOTES.md` — exposition and revision history for the final draft.
- `notes/PROVENANCE.json` — provenance metadata for the proof-audit bundle.

## High-level idea

The Thorp shuffle is rewritten as a walk on the binary hypercube. In a rotating frame, each ordinary shuffle updates one binary coordinate; after \(d\) shuffles every coordinate has been updated once, which is called a sweep.

The proof then combines five ingredients:

1. **Collision smoothing.** Track a small ordered set of marked cards through an independent sweep. Collisions in the routing network yield quantitative contraction for tuple kernels.
2. **Entropy dissipation.** A sweep reveals random partners between cards. Relative entropy decreases at a rate controlled by a scalar defect and a collision term. A near-critical moment interpolation sharpens the entropy clock.
3. **Representation amplification.** Bounds for ordered \(h\)-tuples are transferred to Fourier levels of the full permutation representation of \(S_n\).
4. **Positive truncation.** A high-probability regular component is extracted from a block without replacing the actual shuffle law by a conditional kernel.
5. **Asymmetric three-block finish.** Two entropy-regularized blocks control high levels, while a shorter third block supplies the remaining low-level contraction. Optimizing their lengths yields the leading constant \(5.1536951454\ldots\).

The paper is written so that the physical shuffle-to-hypercube correspondence, entropy quantities, representation-theoretic notation, and block construction can be followed from first principles.

## Building the paper

A standard TeX installation with `pdflatex` is sufficient:

```bash
cd paper
pdflatex main.tex
pdflatex main.tex
pdflatex main.tex
```

No external figures, bibliography database, or data files are required. The repository's GitHub Actions workflow also performs this build and uploads `paper/main.pdf` as an artifact.

## Reproducing the audit checks

The verifier uses only the Python standard library plus NumPy for the floating-point regression checks:

```bash
python3 checks/verify_audit.py
```

The exact combinatorial portions use integer/rational arithmetic. The floating-point checks are diagnostics and are recorded separately from the exact checks in the output.

## Status

This repository records a research manuscript and its computational audit material. The theorem is a new research claim and should be independently checked in the usual way before relying on it as an established published result.

## Citation

If you use this work, please cite the manuscript by Yunjiang Jiang, *A subquadratic upper bound for the mixing time of the Thorp shuffle* (2026). An arXiv identifier can be added here after posting.
