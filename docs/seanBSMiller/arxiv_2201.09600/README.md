# Miller PRC 106 024001 (2022) — arxiv:2201.09600 TeX source

Local cache of the arxiv source tarball, kept for verbatim formula extraction.
The published-PDF OCR is unreliable: e.g. Eq. (D4) phase was rendered as
`(-1)^(J+j+l+1/2)` instead of the correct `(-1)^(J+j)`, and the
`sqrt((2 Sigma + 1))` factor was silently dropped. These errors broke our
implementation until the arxiv source was consulted (commit `a4417a9`).

## Contents

- `source.tar.gz` — full upstream tarball from `https://arxiv.org/e-print/2201.09600`
- `article.tex` — main paper source
- `preamble.tex` — `\input{}`-included preamble
- `bibliography.bib`, `article.bbl` — bibliography
- (Figure PDFs are in the tarball, not copied loose to keep this dir lean.)

## Equations of record (Appendix D, near lines 1345-1410 of `article.tex`)

- **Eq. (D2)**: M-matrix partial-wave summation in channel-spin (Sigma) basis, with
  `i^(l'-l)`, channel-spin Clebsch-Gordans, and `Y_{l'}^{m_Sigma - m_Sigma'}(theta, 0)`.
- **Eq. (D3)**: U -> S relation. `S^J = delta - 2 pi i q m_N i^(l'-l) U^J`. The
  `i^(l'-l)` factors in D2 and D3 cancel when U is plugged in directly, so the
  net geometric factor has NO `i^(l'-l)` phase.
- **Eq. (D4)**: jj -> channel-spin recoupling.
  `U^J_{l' Sigma', l Sigma} = sum_{j', j} sqrt(j-hat' Sigma-hat') (-1)^(J+j') {l' 1/2 j'; J_d J Sigma'} * sqrt(j-hat Sigma-hat) (-1)^(J+j) {l 1/2 j; J_d J Sigma} U^J_{l' j', l j}`,
  with `j-hat = 2 j + 1`.
- **Eq. (D5)**: WP -> plane-wave U conversion via `|f-bar(q)|^2 / (q^2 D-bar_j)`.

These are implemented in `examples/pw_amplitudes.py` (commit `a4417a9`).

## How this was retrieved

Subagent run on 2026-05-10 fetched `https://arxiv.org/e-print/2201.09600` and
extracted the source. Direct alternative: `arxiv.org/abs/2201.09600` -> "Other
formats" -> "Download source".
