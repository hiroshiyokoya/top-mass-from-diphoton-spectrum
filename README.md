# top-mass-from-diphoton-spectrum

Calculation of the $gg\to\gamma\gamma$ diphoton mass spectrum near the top-pair threshold,
$m_{\gamma\gamma}\simeq 2m_t$, at hadron colliders, and its use to determine the top-quark mass.

In $gg\to\gamma\gamma$ the top quark enters only through the one-loop box. Its absorptive part opens
at $m_{\gamma\gamma}=2m_t$ and produces a dip–bump structure in $d\sigma/dm_{\gamma\gamma}$. Near
threshold, the top-loop amplitude is matched to the NRQCD Green function of the $t\bar t$
system,

```math
M_t^{\rm match} = M_t^{\rm 1loop} + \mathcal{B}_t \left[ G(\vec 0; E+i\Gamma_t) - G^{(0)}(\vec 0; E) \right] ,
\qquad E = m_{\gamma\gamma} - 2m_t ,
```

where $G^{(0)}$ is the Green function without QCD interaction and width.
This resums the Coulomb (bound-state) effects. The location and shape of the structure depend on $m_t$ and $\Gamma_t$.

- S. Kawabata and H. Yokoya, *Top-quark mass from the diphoton mass spectrum*,
  [arXiv:1607.00990](https://arxiv.org/abs/1607.00990), Eur. Phys. J. C 77 (2017) 323. This
  repository implements that calculation.
- L. Chen, G. Heinrich, S. Jahn, S. P. Jones, M. Kerner, J. Schlenk and H. Yokoya,
  *Photon pair production in gluon fusion: Top quark effects at NLO with threshold matching*,
  [arXiv:1911.09314](https://arxiv.org/abs/1911.09314), JHEP 04 (2020) 115. This is the NLO
  follow-up and is not implemented here.

The repository has two parts:

- **Signal calculation (main), in Fortran.** It computes $d\sigma/dm_{\gamma\gamma}$ with the
  Green-function matching, produces mass templates, and reproduces the figures of
  arXiv:1607.00990.
- **Mass fit (sub), in Python** (`fit/`). Pseudo-experiments that fit templates plus a smooth
  background to the diphoton spectrum. This part comes from the 2018 code, previously the
  repository `pyTMDP`. See [fit/README.md](fit/README.md).

> **Status (2026-10).**
>
> - The **signal of arXiv:1607.00990 is reproduced**: all 19 curves of its Figs. 1 and 4 agree
>   with this repository to within 0.31% (see
>   [Reproducing arXiv:1607.00990](#reproducing-arxiv160700990)).
> - The fit part works with current ROOT. P1, P3, P4 and P6 of [docs/REVIEW.md](docs/REVIEW.md)
>   are fixed; the remaining items there are open.

## Repository layout

| Path | Content |
|---|---|
| `fortran/lib/amp/` | $gg\to\gamma\gamma$ helicity amplitudes: light-quark loop, top loop, threshold part $\mathcal{A}+\mathcal{B}G$ |
| `fortran/lib/green/` | $t\bar t$ Green function: Schrödinger equation with the LO/NLO QCD potential (`GrnMSBNLO.f`), analytic forms (`GrnMSBLO.f`) |
| `fortran/lib/qcd/` | libQCD: running $\alpha_s$ and quark masses |
| `fortran/lib/include/`, `fortran/lib/bases_stub.f` | shared include file; no-op BASES routines |
| `fortran/extern/` | third-party code (VEGAS, QCD-PEGASUS, HPLOG, DDILOG, CHAPLIN header); see [THIRD_PARTY.md](THIRD_PARTY.md) |
| `fortran/src/mktemplate.f` | current driver: $d\sigma/dm_{\gamma\gamma}$ for one $(m_t,\Gamma_t)$, parameters via namelist |
| `fortran/legacy/` | the 2016 drivers of arXiv:1607.00990 (`MKD_gg2aa.f`, `Sig_gg2aa.f`, …), kept as they were apart from the P3 fix |
| `fortran/Makefile` | `make -C fortran` builds `mktemplate`; `make -C fortran legacy` builds the 2016 drivers |
| `scripts/make_templates.py` | runs `mktemplate` over a $(m_t,\Gamma_t)$ grid in parallel |
| `scripts/reproduce_1607_00990.py` | computes all curves of Figs. 1 and 4 of the paper and compares them with `reference/` |
| `scripts/tools/digitize_1607_00990.py` | extracts the curves from the figure PDFs of the arXiv source |
| `config/templates/` | template sets: the paper's setup (`paper1607_*`), the 2018 sets used by `fit/`, a quick test set |
| `config/repro/1607.00990.yml` | setup and curves for reproducing the paper |
| `reference/1607.00990/` | curves extracted from the paper's vector figures (CSV) |
| `docs/` | code review ([REVIEW.md](docs/REVIEW.md)) and reproduction results ([repro-1607.00990/](docs/repro-1607.00990/summary.md)) |
| `fit/` | mass fit (sub): `TMDP.py`, `ScanMass.py`, `ScanWidth.py`, `Scan2D.py`, fit inputs in `fit/config/` |
| `tests/fortran/`, `tests/fit/` | spot checks against the paper; end-to-end smoke test and config consistency of the fit |
| `docker/Dockerfile` | toolchain image `tmdp`: gfortran, LHAPDF 6.5.5 + CT14 sets, CHAPLIN 1.2, ROOT 6.34, poppler |

## Quick start (Docker)

Everything runs in a container, with the repository mounted at `/work`.

```bash
docker build -t tmdp:dev docker                      # once
docker run --rm -v "$PWD:/work" tmdp:dev make -C fortran
docker run --rm -v "$PWD:/work" tmdp:dev python3 -m pytest -q tests
```

On Windows Git Bash, prefix the commands with `MSYS_NO_PATHCONV=1` and use an absolute path such as
`-v "D:/Physics/top-mass-from-diphoton-spectrum:/work"`. For an interactive shell, use
`docker run --rm -it -v "$PWD:/work" tmdp:dev`.

## Signal calculation

`fortran/build/mktemplate.exe` computes $d\sigma(pp\to gg\to\gamma\gamma)/dm_{\gamma\gamma}$ [fb/GeV]
on an $m_{\gamma\gamma}$ grid for one $(m_t,\Gamma_t)$. The amplitude includes:

- the five light-quark loops,
- the top loop,
- the Green-function correction $\mathcal{B}(G-G_0)$.

The $p_T$/$\eta$ integral is done with VEGAS. Parameters are read from a namelist on standard input;
see the header of `fortran/src/mktemplate.f`.

| Namelist | Meaning | Default |
|---|---|---|
| `RS` | $\sqrt s$ [GeV] | 100000 |
| `PDFSET` | LHAPDF set | `CT14lo` |
| `MT`, `GT` | $m_t$, $\Gamma_t$ [GeV] | 173, 1.5 |
| `MUG` | scale $\mu$ of $\alpha_s$ in the Green function [GeV] | 40 |
| `IORDG` | QCD potential in the Green function: 0 LO, 1 NLO | 1 |
| `MODE` | top amplitude: 3 matched (one-loop + $\mathcal{B}(G-G_0)$), 0 one-loop only | 3 |
| `NQCDIN` | order of the $\alpha_s$ running (libQCD), from $\alpha_s(M_Z)=0.1185$ | 0 (LO) |
| `ETMAX`, `PTMIN`, `PTRATIO` | $\lvert\eta_\gamma\rvert<$ `ETMAX`, $p_T^\gamma>$ `PTMIN`, $p_T^\gamma>$ `PTRATIO` $\cdot m_{\gamma\gamma}$ (0 = off) | 2.5, 40, 0.4 |
| `MAAMIN`, `MAAMAX`, `DMAA` | $m_{\gamma\gamma}$ grid [GeV] | 300, 400, 0.1 |
| `NCALL`, `ITMX` | VEGAS calls and iterations per point | 50000, 6 |
| `LEGACYCUT` | 1 = reproduce the 2016 behaviour outside the cuts ([REVIEW](docs/REVIEW.md) P3) | 0 |
| `OUTFILE` | output: `m_aa  dsigma/dm_aa  error` per line | `Tab_<MT>_<GT>.dat` |

Fixed in the code: $\alpha=1/128$ and $\mu_R=\mu_F=m_{\gamma\gamma}$.

### Templates

```bash
docker run --rm -v "$PWD:/work" tmdp:dev \
    python3 scripts/make_templates.py config/templates/paper1607_LHC13.yml -j 16
```

- A template set (`config/templates/*.yml`) fixes the physics settings and the $(m_t,\Gamma_t)$
  points. Points are given as `masses` × `widths`, as a list of such `grids`, or as explicit
  `points`.
- Each $(m_t,\Gamma_t)$ point is written to `outdir` as `Tab_<mt>_<Gt>.dat`, together with a
  `manifest.json` (configuration and git commit).
- Each template is split into chunks of `--chunk` $m_{\gamma\gamma}$ points (default 25). The
  chunks of all templates run as separate `mktemplate` processes, `-j` at a time, and are joined
  afterwards.
  - An interrupted run resumes from the finished chunks.
  - Each chunk first skips the random numbers that the preceding $m_{\gamma\gamma}$ points would
    have used (`NSKIP` in `mktemplate.f`), so the chunks do not repeat each other's VEGAS random
    sequence.
  - When VEGAS uses exactly `NCALL` calls per iteration and runs all `ITMX` iterations, the result
    equals that of one run per template. For example, with VEGAS 5000 × 3 the mean relative
    difference was $3\times10^{-9}$. Otherwise the sequences are only independent, and the
    numbers differ within the integration errors.
  - With `--chunk 1001` (no split) the output is byte-identical to one run per template.
- Existing templates are skipped unless `--force` is given; `--check` only lists the points.
- Cost: with VEGAS 50 000 × 6 at $\sqrt s=13$ TeV, one $m_{\gamma\gamma}$ point takes about 13–15 s
  on one core (measured with 19–20 parallel processes). A template of 1001 points is therefore
  about 4 core-hours. The CPU time grows roughly linearly with `ncall` × `itmx`.

| Template set | Setup | Templates |
|---|---|---|
| `paper1607_LHC13`, `paper1607_FCC100` | arXiv:1607.00990: CT14nlo, $\Gamma_t=1.498$ GeV, LHC $p_T>40$ GeV / FCC $p_T>0.4m_{\gamma\gamma}$; $m_t$ = 165–181 GeV | 33 each |
| `LHC13T`, `LHC13L`, `HELHC27T`, `FCC100` | the 2018 sets read by `fit/config/`: mass scans at $\Gamma_t=1.5$ GeV plus width and 2D scans, CT14lo assumed; T = tight, L = loose $p_T$ cut (inferred) | 50, 33, 185, 379 |
| `test` | two masses, low statistics, for `tests/` | 2 |

## Reproducing arXiv:1607.00990

```bash
docker run --rm -v "$PWD:/work" tmdp:dev make -C fortran
docker run --rm -v "$PWD:/work" tmdp:dev \
    python3 scripts/reproduce_1607_00990.py config/repro/1607.00990.yml -j 20
```

This computes every curve of Figs. 1 and 4 of the paper with the paper's setup: 19 curves, about
6000 $m_{\gamma\gamma}$ points, roughly 40 min on 20 cores.

- **Outputs.** In `results/1607.00990/`:
  - one `.dat` file per curve
  - plots of our curves over the paper's (`fig1L.png`, `fig1R.png`, `fig4L.png`, `fig4R.png`)
  - `summary.md`
- **Comparison.** Each curve is compared with the same curve extracted from the paper's vector
  figures (`reference/1607.00990/`, made by `scripts/tools/digitize_1607_00990.py` from the
  arXiv source). The script exits with status 0 when every curve agrees within 1%.
- **Spot checks.** `tests/fortran/test_repro_1607_00990.py` checks six points (dip, bump, LO Green
  function, FCC) to within 0.5% in about 30 s.

Setup (`config/repro/1607.00990.yml`):

- **From the paper.** CT14NLO; $\mu_R=\mu_F=m_{\gamma\gamma}$; $|\eta_\gamma|<2.5$.
  - Photon cut: $p_T^\gamma>40$ GeV at the LHC, and additionally $p_T^\gamma>0.4m_{\gamma\gamma}$
    at the FCC.
  - Top quark: $m_t=173$ GeV, $\Gamma_t=1.498$ GeV.
  - Green function: NLO (Fig. 1 left: LO), with $\mu=40$ GeV (Fig. 1: 20–160 GeV).
- **Not in the paper, but needed to match it.** $\alpha_s$ is run at LO from
  $\alpha_s(M_Z)=0.1185$ (`NQCDIN=0`). With NLO running the cross section is 1.4% lower everywhere.

**Result (2026-10-08): all 19 curves agree with the paper.** The largest deviation is 0.31%,
and the mean deviation of every curve is below 0.02%. The dip and bump positions agree to within
one grid step (0.1 GeV for Fig. 1, 0.25 GeV for Fig. 4). The full table and the plots are in
[docs/repro-1607.00990/](docs/repro-1607.00990/summary.md).

| Figure | Curves | Max. abs. deviation |
|---|---|---|
| Fig. 1 right (LHC, $G_{\rm NLO}$, $\mu$ = 20–160 GeV, one-loop) | 5 | 0.17% |
| Fig. 1 left (LHC, $G_{\rm LO}$, $\mu$ = 20–160 GeV, one-loop) | 5 | 0.21% |
| Fig. 4 left (LHC, $m_t$ = 167–179 GeV) | 5 | 0.12% |
| Fig. 4 right (FCC, $p_T>0.4m_{\gamma\gamma}$, $m_t$ = 167–179 GeV) | 5 | 0.31% |

![Fig. 1 right reproduced](docs/repro-1607.00990/fig1R.png)

## The 2016 drivers (`fortran/legacy/`)

These are the programs used for arXiv:1607.00990: `MKD_gg2aa.f` (templates), `Sig_gg2aa.f`,
`DSig_*.f`, `Scl_gg2aa.f`, `Argand_gg2aa.f`, and others. They are kept for reference, and
`make -C fortran legacy` builds the seven that do not need BASES.

As written in 2016, their integrands returned an undefined value for phase-space points rejected by
the cuts. Depending on the compiler, the result changed by up to a factor of 2
([REVIEW](docs/REVIEW.md) P3). This is fixed by one line in each integrand (`INT2 = 0D0`); nothing
else in these files is changed. The old behaviour can be reproduced exactly with `mktemplate`
(`LEGACYCUT=1`), and the paper's figures correspond to the fixed behaviour.

## Mass fit (`fit/`)

The 2018 pseudo-experiment code: template fit of $m_t$ and $\Gamma_t$ with a smooth background
function. The background samples are inputs (not produced here). See
[fit/README.md](fit/README.md).

## Licence

[MIT](LICENSE), © 2016–2026 Hiroshi Yokoya. You may use, modify and redistribute the code freely,
as long as the copyright notice and the licence text are kept. If you use it in a publication,
please cite arXiv:1607.00990.

The MIT licence does not cover the third-party code in `fortran/extern/`: HPLOG, DDILOG,
QCD-PEGASUS, a VEGAS variant and the CHAPLIN interface header. These are shipped unchanged and
remain under the terms of their authors; see [THIRD_PARTY.md](THIRD_PARTY.md). Their licences have
not been checked yet.

## Author

Hiroshi Yokoya
