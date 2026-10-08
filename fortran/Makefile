# Linux/Docker build of the gg2aa programs (arXiv:1607.00990) and libQCD.
# Run inside the pytmdp Docker image:  make -C fortran [target]
# The original macOS makefiles are kept as */Makefile.legacy.

FC       = gfortran
FFLAGS  ?= -O2 -w -std=legacy -ffixed-line-length-132
BUILD    = build
INCS     = -IQCD -Igg2aa -I$(BUILD)/inc

CHAPLIN_DIR ?= /opt/chaplin
LHAPDF_LIBS  = $(shell lhapdf-config --ldflags) -lstdc++
LIBS         = $(BUILD)/libQCD.a -L$(CHAPLIN_DIR)/lib -lchaplin $(LHAPDF_LIBS)

QCD_OBJS   = $(addprefix $(BUILD)/,QCDlibs.o asrgkt_pegasus.o betafct_pegasus.o)
GG2AA_OBJS = $(addprefix $(BUILD)/,gg2aaQ.o gg2aaT.o gg2aaG.o gg2aaA.o \
               GrnMSBLO.o GrnMSBNLO.o intvegas.o bases_stub.o)

# Sig_tot.f integrates with BASES itself and is not built (BASES is not
# shipped; fortran/stubs only provides no-op set-up/histogram routines).
PROGRAMS = mktemplate MKD_gg2aa Sig_gg2aa Sig_qqb2aa Amp_gg2aa Argand_gg2aa \
           DGrnMS PlotGrnMS
EXES     = $(addprefix $(BUILD)/,$(addsuffix .exe,$(PROGRAMS)))

.PHONY: all clean
.SECONDARY:
all: $(EXES)

# gg2aaT.f / gg2aaA.f include 'CHAPLIN.INC' but the file is shipped as
# chaplin.inc; provide the upper-case name on case-sensitive file systems.
$(BUILD)/inc/CHAPLIN.INC: gg2aa/chaplin.inc
	@mkdir -p $(BUILD)/inc
	cp $< $@

$(BUILD)/libQCD.a: $(QCD_OBJS)
	ar cr $@ $^
	ranlib $@

$(BUILD)/%.o: QCD/%.f | $(BUILD)/inc/CHAPLIN.INC
	$(FC) $(FFLAGS) $(INCS) -c $< -o $@

$(BUILD)/%.o: gg2aa/%.f | $(BUILD)/inc/CHAPLIN.INC
	$(FC) $(FFLAGS) $(INCS) -c $< -o $@

$(BUILD)/%.o: src/%.f | $(BUILD)/inc/CHAPLIN.INC
	$(FC) $(FFLAGS) $(INCS) -c $< -o $@

$(BUILD)/%.o: stubs/%.f
	@mkdir -p $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/%.exe: $(BUILD)/%.o $(GG2AA_OBJS) $(BUILD)/libQCD.a
	$(FC) $(FFLAGS) $< $(GG2AA_OBJS) $(LIBS) -o $@

clean:
	rm -rf $(BUILD)
