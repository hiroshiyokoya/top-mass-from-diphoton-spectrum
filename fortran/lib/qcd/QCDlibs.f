C     QCDlibs.f (libQCD): running alpha_s and quark-mass conversions.
C     alpha_s is evolved with QCD-PEGASUS (fortran/extern/pegasus).
C     COMMON blocks (qcdparam.inc):
C       /QCDCON/ PI
C       /QCDIN/  MZ, ASZ = alpha_s(MZ), NQCD
C       /PMASS/  pole masses PMT, PMB, PMC, PMS [GeV]
C       /MASS0/  MSbar masses MT0, MB0, MC0, MS0 [GeV]
C       /QCDLIB/ NC, NF, CA, CF, TR (B0, QTHR are not set here)
C     NQCD is the order of the alpha_s running, passed to PEGASUS as
C     NAORD (0: LO, solved exactly; 1: NLO, 2, 3: Runge-Kutta).
C     It must be set before the first ASQCD call (QCDLIBS makes one),
C     because PEGASUS is set up only once.
C
C     QCDLIBS: sets MZ = 91.2 GeV, alpha_s(MZ) = 0.1185, n_f = 5, the
C     colour factors and the quark masses: top pole mass 173.2 GeV
C     (-> MSbar MT0); MSbar m_b = 4.2, m_c = 1.275, m_s = 0.115 GeV
C     (-> pole masses).  PMB, PMC and PMT are the flavour thresholds of
C     ASQCD.  PMT is fixed here and does not follow the m_t of a run.
C     The conversions call ASQCD while PMB and PMC are still being set
C     (not yet set), so each uses the thresholds set so far.
      SUBROUTINE QCDLIBS
      IMPLICIT NONE
      INCLUDE 'qcdparam.inc'
      INCLUDE 'qcdfunc.inc'
C...  Constant
      PI  = DACOS(-1D0)
C...  Weak bosons
      MZ   = 91.2D0             ! 
C.....
      NC   = 3D0
      CA   = NC
      CF   = (NC**2-1D0)/(2D0*NC)
      TR   = 0.5D0
C...  QCD coupling
      NF    = 5
      ASZ   = 0.1185D0          ! 
C...  Input Quark Masses
      PMT  = 173.2D0            ! Top     Pole mass
      MT0  = MSBARMASS(PMT)
      MB0  = 4.2D0              ! Bottom  MSbar mass
      PMB  = POLEMASS (MB0)
      MC0  = 1.275D0            ! Charm   MSbar mass
      PMC  = POLEMASS (MC0)
      MS0  = 0.115D0            ! Strange MSbar mass
      PMS  = POLEMASS (MS0)
C.....
      RETURN
      END
C.....
C
C     ASQCD (Q): MSbar alpha_s at the scale Q [GeV], from alpha_s(MZ),
C     evolved by AS_PEGASUS at order NQCD (in a = alpha_s/(4 pi)) in a
C     variable-flavour scheme:
C       Q > PMT        : n_f = 5 from PMT (the n_f = 6 line is commented
C                        out, so this equals 5-flavour running from MZ)
C       PMB < Q <= PMT : n_f = 5 from MZ
C       PMC < Q <= PMB : n_f = 4 from PMB
C       Q <= PMC       : n_f = 3 from PMC
C     At the b and c thresholds alpha_s is matched by the factor
C     (1 - 14/3 a**2); AST6 is computed but not used.  The Q.EQ.Q0
C     branch has no effect (Q0 is never updated and there is no RETURN).
      DOUBLE PRECISION FUNCTION ASQCD (Q)
      IMPLICIT NONE
      DOUBLE PRECISION Q0,AS0,Q,Q1,AS_PEGASUS
      COMMON /QCD0/ Q0,AS0
      DATA Q0 /1.234D5/
      INTEGER IASQCD
      DATA    IASQCD /0/
      DOUBLE PRECISION AST6,AST5,ASB5,ASB4,ASC4,ASC3
      COMMON /PEGASUS/ AST6,AST5,ASB5,ASB4,ASC4,ASC3,IASQCD
      DOUBLE PRECISION CF_PEGASUS,CA_PEGASUS,TR_PEGASUS
      COMMON /COLOUR/  CF_PEGASUS,CA_PEGASUS,TR_PEGASUS
      INTEGER        NAORD,NASTPS
      COMMON /ASPAR/ NAORD,NASTPS
      INCLUDE 'qcdparam.inc'
      IF ( Q.EQ.Q0 ) THEN
         ASQCD = AS0
      ENDIF
      IF ( IASQCD.EQ.0 ) THEN
         CALL INITASQCD
         IASQCD = 1
      ENDIF
      AST5 = AS_PEGASUS(PMT**2, MZ**2,ASZ /(4D0*PI),5)*(4D0*PI)
      AST6 = AST5 * ( 1D0 + 14D0/3D0*AST5**2/(4D0*PI)**2 )
      ASB5 = AS_PEGASUS(PMB**2, MZ**2,ASZ /(4D0*PI),5)*(4D0*PI)
      ASB4 = ASB5 * ( 1D0 - 14D0/3D0*ASB5**2/(4D0*PI)**2 )
      ASC4 = AS_PEGASUS(PMC**2,PMB**2,ASB4/(4D0*PI),4)*(4D0*PI)
      ASC3 = ASC4 * ( 1D0 - 14D0/3D0*ASC4**2/(4D0*PI)**2 )
      Q1 = Q
c     IF ( Q .LT.1D0 ) Q1 = 1D0
      IF ( Q1.GT.PMT ) THEN
c     ASQCD = AS_PEGASUS(Q1**2,PMT**2,AST6/(4D0*PI),6)*(4D0*PI)
         ASQCD = AS_PEGASUS(Q1**2,PMT**2,AST5/(4D0*PI),5)*(4D0*PI)
      ELSEIF ( Q1.GT.PMB ) THEN
         ASQCD = AS_PEGASUS(Q1**2, MZ**2,ASZ /(4D0*PI),5)*(4D0*PI)
      ELSEIF ( Q1.GT.PMC ) THEN
         ASQCD = AS_PEGASUS(Q1**2,PMB**2,ASB4/(4D0*PI),4)*(4D0*PI)
      ELSE
         ASQCD = AS_PEGASUS(Q1**2,PMC**2,ASC3/(4D0*PI),3)*(4D0*PI)
      ENDIF
      AS0 = ASQCD
      RETURN
      END
C.....
C
C     INITASQCD: sets up PEGASUS: NAORD = NQCD, NASTPS = 10 Runge-Kutta
C     steps, colour factors, beta-function coefficients
C     (BETAFCT_PEGASUS); prints the order.
      SUBROUTINE INITASQCD
      IMPLICIT NONE
      DOUBLE PRECISION CF_PEGASUS,CA_PEGASUS,TR_PEGASUS
      COMMON /COLOUR/  CF_PEGASUS,CA_PEGASUS,TR_PEGASUS
      INTEGER        NAORD,NASTPS
      COMMON /ASPAR/ NAORD,NASTPS
      INCLUDE 'qcdparam.inc'
      NAORD  = NQCD
      NASTPS = 10
      CF_PEGASUS = CF
      CA_PEGASUS = CA
      TR_PEGASUS = TR
      CALL BETAFCT_PEGASUS
      PRINT *, "Set QCD-Pegasus with N^nLO expansion with n = ", NQCD
      RETURN
      END
C.....
C
C     POLEMASS (M0): pole mass from the MSbar mass M0 [GeV] at one loop,
C     m_pole = M0 (1 + C_F alpha_s(M0)/pi).
      DOUBLE PRECISION FUNCTION POLEMASS (M0)
      IMPLICIT NONE
      DOUBLE PRECISION M0, ASQCD,ASM,D1
      EXTERNAL ASQCD
      INCLUDE 'qcdparam.inc'
      ASM      = ASQCD(M0)
      D1       = CF
      PoleMass = M0*(1D0+ASM/PI*D1)
      RETURN
      END
C.....
C
C     MSBARMASS (M0): inverse of POLEMASS for a pole mass M0 [GeV],
C     m = M0/(1 + C_F alpha_s(M0)/pi).  The NNLO variant is commented
C     out.
      DOUBLE PRECISION FUNCTION MSBARMASS (M0)
      IMPLICIT NONE
      DOUBLE PRECISION M0, ASQCD,ASM,D1
      EXTERNAL ASQCD
      INCLUDE 'qcdparam.inc'
      ASM       = ASQCD(M0)
      D1        = CF
      MSBARMASS = M0/( 1D0 + ASM/PI*D1 ) ! NLO
c     MSBARMASS = M0/( 1D0 + ASM/PI*D1 + (ASM/PI)**2*10.9 ) ! NNLO for top [HDECAY]
      RETURN
      END
