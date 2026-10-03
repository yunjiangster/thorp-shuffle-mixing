# Expository revision notes

This version preserves the theorem, constants, lemmas, and detailed proofs from the audited arXiv draft, while substantially expanding the exposition for readers new to the Thorp shuffle and to mixing-time/representation-theoretic arguments.

Major changes:

- Added a step-by-step derivation of the hypercube model from the physical two-packet shuffle.
- Explicitly factorized one physical shuffle as a random coordinate-1 matching layer followed by cyclic coordinate rotation.
- Proved the rotating-frame identity showing that coordinate directions cycle and that one sweep equals exactly d ordinary Thorp shuffles.
- Added a complete n=8 example, including the possible motion of a card initially at binary position 101.
- Explained why one card is exactly uniform after one sweep and why full-deck mixing remains difficult because of dependence among cards.
- Added an explicit definition of worst-start total-variation mixing time and of log^+.
- Expanded the representation-theory notation to define partitions, Young diagrams, conjugation, levels, and standard Young tableaux.
- Added a first-time-reader dictionary connecting entropy, marked cards, tuple kernels, representation levels, truncation, and shuffle blocks.
- Rewrote the proof overview around the geometric and global tasks and added a roadmap table.
- Added plain-language introductions to the entropy clock, ordered-tuple amplification, truncation, asymmetric three-block finish, and coefficient optimization.
- Added orientation paragraphs to each technical appendix explaining why the next calculation is needed and how it connects to the main proof.
- Clarified the practical consequence of the two external inputs before stating their formal versions.
- Replaced the phrase “individual layers” in the final finite criterion by “ordinary Thorp shuffles” to keep the physical clock unambiguous.

The mathematical result remains

    t_mix(epsilon) <= 5.153695145414931... d^(5/3)
                      + O_epsilon(d^(4/3) log(d+2)),

for n=2^d, with physical time measured in ordinary Thorp shuffles.

The final PDF has 31 pages and was compiled repeatedly with pdflatex. PDF preflight passed and the rendered pages were visually inspected, including the new hypercube and roadmap sections.
