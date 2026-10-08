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
      DOUBLE PRECISION FUNCTION POLEMASS (M0)
      IMPLICIT NONE
      DOUBLE PRECISION M0,MP,ASQCD,ASM,D1
      EXTERNAL ASQCD
      INCLUDE 'qcdparam.inc'
      ASM      = ASQCD(M0)
      D1       = CF
      PoleMass = M0*(1D0+ASM/PI*D1)
      RETURN
      END
C.....
      DOUBLE PRECISION FUNCTION MSBARMASS (M0)
      IMPLICIT NONE
      DOUBLE PRECISION M0,MP,ASQCD,ASM,D1
      EXTERNAL ASQCD
      INCLUDE 'qcdparam.inc'
      ASM       = ASQCD(M0)
      D1        = CF
      MSBARMASS = M0/( 1D0 + ASM/PI*D1 ) ! NLO
c     MSBARMASS = M0/( 1D0 + ASM/PI*D1 + (ASM/PI)**2*10.9 ) ! NNLO for top [HDECAY]
      RETURN
      END
