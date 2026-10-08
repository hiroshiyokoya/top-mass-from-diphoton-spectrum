# Mass fit (`fit/`)

Pseudo-experiments for the top-quark mass (and width) determination from the diphoton spectrum, as
in Sec. 5 of [arXiv:1607.00990](https://arxiv.org/abs/1607.00990). This is the 2018 Python code,
previously the whole repository `pyTMDP`. It needs signal templates from the Fortran part (see the
[main README](../README.md)) and background samples, which are inputs.

> **Status (2026-10).** P1 (empty pseudo-data with current ROOT), P4 and P6 (template names parsed
> from the full path) of [docs/REVIEW.md](../docs/REVIEW.md) are fixed. The pseudo-data are now
> Poisson numbers around the expected counts per bin; set `seed` in the fit input to make a run
> reproducible.
>
> Still open:
>
> - **P2** The best-fit mass is the template with the smallest $\chi^2$, so it is quantised to the
>   template grid.
> - **P5** The fit function takes the template value at the bin centre instead of the bin
>   integral.
> - **P20** The background shape in the pseudo-data is the histogram of the background sample
>   itself, so the sample must be much larger than one pseudo-dataset. For example, 20 million
>   events give $\chi^2/{\rm ndf}\simeq1$, but 0.2 million give about 2.5.

## How it works

1. **Signal templates.** One file `Tab_<mt>_<Gt>.dat` per $(m_t,\Gamma_t)$ holds
   $d\sigma/dm_{\gamma\gamma}$ for 300–400 GeV in 0.1 GeV steps, as written by
   `scripts/make_templates.py`.
2. **Backgrounds.** Event lists of $m_{\gamma\gamma}$ for direct, one-fragmentation and
   two-fragmentation photon pairs (`Direct.dat`, `OneF.dat`, `TwoF.dat`).
3. **Pseudo-experiments** (`TMDP.py`). The background and the signal template at the true mass
   are integrated over `hbin` bins of 300–400 GeV, scaled to the expected numbers of events, and the
   observed number in each bin is drawn from a Poisson distribution.
   - Each pseudo-dataset is fitted in `[fmin, fmax]` with
     $(1-k_{gg})\,f_{\rm ATLAS}(m;a)+k_{gg}\,f_{\rm template}(m;m_t)$, where
     $f_{\rm ATLAS}\propto(1-(m/\sqrt s)^{1/3})^a$, once for every template.
   - The template with the smallest $\chi^2$ gives the best-fit $m_t$.
   - `ScanMass.py` repeats this `Nloop` times and reports the mean and spread. `ScanWidth.py` and
     `Scan2D.py` do the same for $\Gamma_t$ and for $(m_t,\Gamma_t)$.

## Running

Make the templates, put the background files in the same directory, then run from the repository
root (paths in the fit inputs are relative to it):

```bash
docker run --rm -v "$PWD:/work" tmdp:dev \
    python3 scripts/make_templates.py config/templates/LHC13T.yml -j 16
docker run --rm -v "$PWD:/work" tmdp:dev python3 fit/ScanMass.py fit/config/input_LHC13T.yml 1000
```

## Fit inputs (`fit/config/`)

The eight inputs of 2018. Each reads its templates and backgrounds from `dir`, which is the `outdir`
of exactly one template set in `config/templates/`. `tests/fit/test_config.py` checks that every
template a fit input needs is in that set.

| Fit input | Template set (`dir`) | $\sqrt s$ | Scan |
|---|---|---|---|
| `input_LHC13T.yml`, `inputWidth_LHC13T.yml` | `Template/LHC13T` | 13 TeV | $m_t$, $\Gamma_t$ |
| `input_LHC13L.yml` | `Template/LHC13L` | 13 TeV | $m_t$ |
| `input_HELHC27T.yml`, `input2D_HELHC27T.yml` | `Template/HELHC27T` | 27 TeV | $m_t$, $(m_t,\Gamma_t)$ |
| `input_FCC100.yml`, `inputWidth_FCC100.yml`, `input2D_FCC100.yml` | `Template/FCC100` | 100 TeV | $m_t$, $\Gamma_t$, $(m_t,\Gamma_t)$ |

Keys, as read by `TMDP.set_init`:

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
| `seed` | optional seed of the pseudo-data random numbers |

The units of `sig_*` are inferred from the `Nevnt` formula and have not been checked against the
original setup.

## Background data

`TMDP.read_events_from_file` expects plain text with one $m_{\gamma\gamma}$ value [GeV] per event,
separated by whitespace, in 300–400 GeV. The paper generated them with Diphox at LO
($q\bar q\to\gamma\gamma$, one- and two-fragmentation; Sec. 5 of arXiv:1607.00990). The 2018 files
are not available.
