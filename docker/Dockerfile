# pyTMDP toolchain: ROOT (PyROOT) + gfortran + LHAPDF6 + CHAPLIN + PDF sets.
# The repository itself is mounted at /work at run time (see README).
FROM rootproject/root:6.34.00-ubuntu24.04

ARG LHAPDF_VERSION=6.5.5
ARG CHAPLIN_VERSION=1.2
ARG PDF_SETS="CT14lo CT14nlo"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        gfortran g++ make wget ca-certificates parallel \
        python3-numpy python3-matplotlib python3-yaml python3-pytest \
    && rm -rf /var/lib/apt/lists/*

# LHAPDF6 (C++ library + LHAGlue Fortran interface; Python bindings not needed)
RUN wget -q -O /tmp/lhapdf.tar.gz "https://lhapdf.hepforge.org/downloads/?f=LHAPDF-${LHAPDF_VERSION}.tar.gz" \
    && tar -xzf /tmp/lhapdf.tar.gz -C /tmp \
    && cd /tmp/LHAPDF-${LHAPDF_VERSION} \
    && ./configure --prefix=/opt/lhapdf --disable-python \
    && make -j"$(nproc)" && make install \
    && rm -rf /tmp/lhapdf.tar.gz /tmp/LHAPDF-${LHAPDF_VERSION}
ENV PATH=/opt/lhapdf/bin:$PATH \
    LD_LIBRARY_PATH=/opt/lhapdf/lib \
    LHAPDF_DATA_PATH=/opt/lhapdf/share/LHAPDF

# PDF sets used by the template generator
RUN for s in ${PDF_SETS}; do \
        wget -q -O - "https://lhapdfsets.web.cern.ch/lhapdfsets/current/${s}.tar.gz" \
        | tar -xz -C /opt/lhapdf/share/LHAPDF; \
    done

# CHAPLIN (complex harmonic polylogarithms, HPL2 used by gg2aaT/gg2aaA)
RUN wget -q -O /tmp/chaplin.tar "https://chaplin.hepforge.org/code/chaplin-${CHAPLIN_VERSION}.tar" \
    && tar -xf /tmp/chaplin.tar -C /tmp \
    && cd /tmp/chaplin-${CHAPLIN_VERSION} \
    && ./configure --prefix=/opt/chaplin \
    && make -j"$(nproc)" && make install \
    && rm -rf /tmp/chaplin.tar /tmp/chaplin-${CHAPLIN_VERSION}
ENV LD_LIBRARY_PATH=/opt/chaplin/lib:/opt/lhapdf/lib \
    CHAPLIN_DIR=/opt/chaplin \
    LHAPDF_DIR=/opt/lhapdf

# poppler-utils: pdftocairo for scripts/tools/digitize_1607_00990.py
RUN apt-get update && apt-get install -y --no-install-recommends poppler-utils \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /work
CMD ["bash"]
