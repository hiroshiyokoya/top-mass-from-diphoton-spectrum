# pyTMDP

**T**op-quark **M**ass from the **D**i**P**hoton mass spectrum.

Pseudo-experiment study of how well the top-quark mass $m_t$ (and width $\Gamma_t$) can be
determined from the shape of the diphoton invariant-mass spectrum $m_{\gamma\gamma}$ near the
top-pair threshold $m_{\gamma\gamma}\simeq 2m_t$ at hadron colliders.

In $gg\to\gamma\gamma$ the top quark enters only through the one-loop box. Its absorptive part opens
at $m_{\gamma\gamma}=2m_t$ and produces a dip–bump structure in $d\sigma/dm_{\gamma\gamma}$. Near
threshold, the amplitude is matched to the NRQCD Green function of the $t\bar t$ system,
$\mathcal{M}^{\rm 1loop}+\mathcal{B}\,[G(E+i\Gamma_t)-G_0(E)]$, which resums the Coulomb (bound-state)
effects. The location and shape of the structure depend on $m_t$ and $\Gamma_t$.

- S. Kawabata and H. Yokoya, *Top-quark mass from the diphoton mass spectrum*,
  [arXiv:1607.00990](https://arxiv.org/abs/1607.00990), Eur. Phys. J. C 77 (2017) 323. This is the
  LO + Green-function calculation that `fortran/gg2aa` implements.
- L. Chen, G. Heinrich, S. Jahn, S. P. Jones, M. Kerner, J. Schlenk and H. Yokoya,
  *Photon pair production in gluon fusion: Top quark effects at NLO with threshold matching*,
  [arXiv:1911.09314](https://arxiv.org/abs/1911.09314), JHEP 04 (2020) 115. This is the NLO
  follow-up and is not implemented here.

> **Status (2026-10).** The 2018 Python code runs with Python 3.12 / ROOT 6.34 in the Docker image,
> but has known problems that are documented and not yet fixed (see [docs/REVIEW.md](docs/REVIEW.md)):
>
> - **P1** With current ROOT, `TH1::FillRandom` silently produces *empty* pseudo-data in
>   `TMDP.genEvents`, so fits run on empty histograms.
> - **P3** The legacy Fortran integrands return an undefined value outside the cuts. Depending
>   on the compiler optimisation, $d\sigma/dm_{\gamma\gamma}$ changes by a factor of about 2.
>   `fortran/src/mktemplate.f`, which generates the templates, is fixed.

## How it works

1. **Templates (signal).** For each $(m_t,\Gamma_t)$, `fortran/build/mktemplate.exe` computes
   $d\sigma(pp\to gg\to\gamma\gamma)/dm_{\gamma\gamma}$ [fb/GeV] for 300–400 GeV in 0.1 GeV steps.
   It includes the light-quark and top loops and the Green-function matching, and writes
   `Tab_<mt>_<Gt>.dat`.
2. **Backgrounds.** These are $m_{\gamma\gamma}$ event lists for direct, one-fragmentation and
   two-fragmentation photon pairs (`Direct.dat`, `OneF.dat`, `TwoF.dat`). They are **not** in this
   repository; see [Background data](#background-data).
3. **Pseudo-experiments** (`TMDP.py`). Events are generated from the background and from the
   signal template at the true mass, and binned in `hbin` bins over 300–400 GeV.
   - Each pseudo-dataset is fitted in `[fmin, fmax]` with
     $(1-k_{gg})\,f_{\rm ATLAS}(m;a)+k_{gg}\,f_{\rm template}(m;m_t)$, where
     $f_{\rm ATLAS}\propto(1-(m/\sqrt s)^{1/3})^a$, once for every template.
   - The template with the smallest $\chi^2$ gives the best-fit $m_t$.
   - `ScanMass.py` repeats this `Nloop` times and reports the mean and spread.

## Repository layout

| Path | Content |
|---|---|
| `TMDP.py` | Classes `TMDP` (inputs, pseudo-data, fit functions) and `GG2AA` (one template) |
| `ScanMass.py`, `ScanWidth.py`, `Scan2D.py` | Scans in $m_t$, in $\Gamma_t$, and in $(m_t,\Gamma_t)$ |
| `input*.yml` | Fit inputs for LHC 13 TeV, HE-LHC 27 TeV and FCC 100 TeV (2018) |
| `fortran/gg2aa/` | 2016 Fortran of arXiv:1607.00990 (amplitudes, Green functions, drivers), unchanged |
| `fortran/QCD/` | libQCD: running $\alpha_s$ (QCD-PEGASUS), unchanged |
| `fortran/src/mktemplate.f` | Template generator, one $(m_t,\Gamma_t)$ per run, parameters via namelist |
| `fortran/stubs/` | No-op BASES routines (BASES is not needed for the VEGAS-based programs) |
| `fortran/Makefile` | Linux build (the original macOS makefiles are kept as `*/Makefile.legacy`) |
| `scripts/make_templates.py` | Runs `mktemplate.exe` over a $(m_t,\Gamma_t)$ grid in parallel |
| `config/templates_*.yml` | Template grids matching `input_*.yml`, plus a quick test grid |
| `docker/Dockerfile` | Toolchain image: ROOT 6.34, gfortran, LHAPDF 6.5.5 + CT14 sets, CHAPLIN 1.2 |
| `tests/` | End-to-end smoke test (Fortran → templates → fit) |
| `docs/REVIEW.md` | Code review (2026-10) and proposed fixes |
| `THIRD_PARTY.md` | Third-party Fortran shipped in `fortran/` |

## Quick start (Docker)

Everything runs in a container, with the repository mounted at `/work`.

```bash
docker build -t pytmdp:dev docker                    # once
docker run --rm -v "$PWD:/work" pytmdp:dev make -C fortran
docker run --rm -v "$PWD:/work" pytmdp:dev python3 -m pytest -q tests
```

On Windows Git Bash, prefix the commands with `MSYS_NO_PATHCONV=1` and use an absolute path such as
`-v "D:/Physics/pyTMDP:/work"`. For an interactive shell, use
`docker run --rm -it -v "$PWD:/work" pytmdp:dev`.

### Generate templates

```bash
docker run --rm -v "$PWD:/work" pytmdp:dev \
    python3 scripts/make_templates.py config/templates_LHC13T.yml -j 16
```

- Templates are written to the `outdir` of the config (for example `Template/LHC13T/`), together
  with a `manifest.json` that records the configuration and the git commit.
- Existing files are skipped unless you pass `--force`.
- One production template (1001 points, VEGAS 50 000 calls × 6 iterations) is estimated to take about 80 min on
  one core (scaled from a timed low-statistics run). `config/templates_LHC13T.yml` has 33 masses, i.e. about 2.5 h on 16–20 cores.

Config keys (all optional except `outdir`, `masses`, `widths`; defaults are those of
`MKD_gg2aa.f`):

| Key | Meaning | Default |
|---|---|---|
| `rs` | $\sqrt s$ [GeV] | 100000 |
| `pdfset` | LHAPDF set (must be installed in the image) | `CT14lo` |
| `masses`, `widths` | list, or `{start, stop, step}` [GeV] | — |
| `mu_green` | scale of $\alpha_s$ in the Green function [GeV] | 40 |
| `etamax`, `ptmin` | photon cuts; in addition $p_T>0.4\,m_{\gamma\gamma}$ is hard-coded | 2.5, 40 |
| `maa_min`, `maa_max`, `maa_step` | $m_{\gamma\gamma}$ grid [GeV]. `TMDP.py` requires 300, 400, 0.1 | 300, 400, 0.1 |
| `ncall`, `itmx` | VEGAS calls per point and iterations | 50000, 6 |

The physics settings are those of `MKD_gg2aa.f`: CT14lo, $\alpha=1/128$, $\mu_R=\mu_F=m_{\gamma\gamma}$,
LO running of $\alpha_s$, and an NLO Coulomb potential in the Green function.

### Run a mass scan

Put the background files in the template directory, then run:

```bash
docker run --rm -v "$PWD:/work" pytmdp:dev python3 ScanMass.py input_LHC13T.yml 1000
```

Fit inputs (`input_*.yml`, as read by `TMDP.set_init`):

| Key | Meaning |
|---|---|
| `dir` | directory with templates and background files |
| `rs`, `lum`, `corr` | $\sqrt s$ [GeV], luminosity [fb⁻¹], overall correction factor |
| `kgg` | $gg\to\gamma\gamma$ (signal) fraction of all events |
| `sig_dir`, `sig_one`, `sig_two` | background cross sections in 300–400 GeV [pb] |
| `Nevnt` | total number of events; if omitted, $\sigma_{\rm bg}\cdot L\cdot 10^3\cdot{\rm corr}/(1-k_{gg})$ |
| `hbin`, `fmin`, `fmax`, `fitopt` | histogram bins in 300–400 GeV, fit range [GeV], ROOT fit options |
| `files_dir`, `files_one`, `files_two` | background event files |
| `files_sig` | template used to generate the pseudo-data (true mass) |
| `files_template` | templates fitted to each pseudo-dataset |

The units of `sig_*` are inferred from the `Nevnt` formula and have not been checked against the
original setup.

## Background data

`TMDP.read_events_from_file` expects plain text with one $m_{\gamma\gamma}$ value [GeV] per event,
separated by whitespace, in 300–400 GeV. The 2018 files are not available, and the generator used
to make them (for example DIPHOX) has not been confirmed. Producing them is outside the scope of the
current milestone.

## Known issues and plans

See [docs/REVIEW.md](docs/REVIEW.md) for the full list and proposed fixes. Fixes are made after
review, in separate issues and PRs. Progress is tracked in issue #7.

## Third-party code

`fortran/` ships HPLOG, DDILOG, QCD-PEGASUS, a VEGAS variant and the CHAPLIN interface header,
unchanged. See [THIRD_PARTY.md](THIRD_PARTY.md). Their licences have not been checked yet.

## Author

Hiroshi Yokoya
