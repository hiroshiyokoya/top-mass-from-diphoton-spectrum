# Third-party code

`fortran/` contains code written by others, shipped unchanged because the
gg2aa programs were built with it in 2016. Their licences have **not** been
checked yet; please consult the original sources before reusing them outside
this repository. They are not covered by the MIT licence of this repository
(`LICENSE`).

| File | Origin | Reference | Licence |
|---|---|---|---|
| `fortran/gg2aa/hplog.f` | HPLOG v1.1 (2004), T. Gehrmann and E. Remiddi | arXiv:hep-ph/0107173, Comput. Phys. Commun. 141 (2001) | to be checked (not linked into the current programs) |
| `fortran/gg2aa/ddilog.f` | `DDILOG`, as in CERNLIB (C332) | — | to be checked (not linked) |
| `fortran/QCD/asrgkt_pegasus.f`, `fortran/QCD/betafct_pegasus.f` | QCD-PEGASUS, A. Vogt | arXiv:hep-ph/0408244, Comput. Phys. Commun. 170 (2005) | to be checked |
| `fortran/gg2aa/intvegas.f` | VEGAS (G. P. Lepage) variant labelled "integration routine for MAX"; random numbers after Knuth | — | origin and licence to be checked |
| `fortran/gg2aa/chaplin.inc` | Interface declarations of CHAPLIN | arXiv:1106.5739, Comput. Phys. Commun. 185 (2014) | to be checked |

Downloaded at image build time and **not** shipped: ROOT (base image
`rootproject/root`), LHAPDF 6, the CT14 PDF sets and CHAPLIN 1.2
(see `docker/Dockerfile`).
