# Build of the gg -> gamma gamma threshold code (arXiv:1607.00990).
# Run inside the tmdp Docker image:  make -C fortran [target]
#
#   lib/      physics library: amplitudes (amp), Green functions (green),
#             running alpha_s (qcd), shared include files, BASES stubs
#   extern/   third-party code (THIRD_PARTY.md), shipped unchanged
#   src/      current drivers (mktemplate)
#   legacy/   the 2016 drivers, kept unchanged; built by `make legacy`

FC       = gfortran
FFLAGS  ?= -O2 -w -std=legacy -ffixed-line-length-132
BUILD    = build
INCS     = -Ilib/qcd -Ilib/include -Iextern -I$(BUILD)/inc

CHAPLIN_DIR ?= /opt/chaplin
LHAPDF_LIBS  = $(shell lhapdf-config --ldflags) -lstdc++
LIBS         = $(BUILD)/libtmdp.a -L$(CHAPLIN_DIR)/lib -lchaplin $(LHAPDF_LIBS)

LIB_SRCS = $(wildcard lib/amp/*.f lib/green/*.f lib/qcd/*.f) lib/bases_stub.f \
           extern/intvegas.f $(wildcard extern/pegasus/*.f)
LIB_OBJS = $(addprefix $(BUILD)/,$(notdir $(LIB_SRCS:.f=.o)))

PROGRAMS = mktemplate
# Sig_tot.f integrates with BASES itself and is not built (BASES is not
# shipped; lib/bases_stub.f only provides no-op set-up/histogram routines).
LEGACY   = MKD_gg2aa Sig_gg2aa Sig_qqb2aa Amp_gg2aa Argand_gg2aa DGrnMS PlotGrnMS

vpath %.f src legacy lib lib/amp lib/green lib/qcd extern extern/pegasus

.PHONY: all legacy clean
.SECONDARY:
all: $(addprefix $(BUILD)/,$(addsuffix .exe,$(PROGRAMS)))
legacy: $(addprefix $(BUILD)/,$(addsuffix .exe,$(LEGACY)))

# lib/amp/gg2aaT.f and gg2aaA.f include 'CHAPLIN.INC' but the file is
# shipped as extern/chaplin.inc; provide the upper-case name on
# case-sensitive file systems.
$(BUILD)/inc/CHAPLIN.INC: extern/chaplin.inc
	@mkdir -p $(BUILD)/inc
	cp $< $@

$(BUILD)/%.o: %.f | $(BUILD)/inc/CHAPLIN.INC
	$(FC) $(FFLAGS) $(INCS) -c $< -o $@

$(BUILD)/libtmdp.a: $(LIB_OBJS)
	ar cr $@ $^
	ranlib $@

$(BUILD)/%.exe: $(BUILD)/%.o $(BUILD)/libtmdp.a
	$(FC) $(FFLAGS) $< $(LIBS) -o $@

clean:
	rm -rf $(BUILD)
