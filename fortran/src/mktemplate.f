C     mktemplate: d sigma / d M_aa template for pp -> gamma gamma via gg
C     near the top-loop threshold (arXiv:1607.00990), one (m_t, Gamma_t)
C     point per run.  Physics and integrand are taken over unchanged
C     from gg2aa/MKD_gg2aa.f; only the run parameters are read from a
C     namelist on standard input instead of being hard-coded.
C
C     Input (namelist /TMPL/, all optional, defaults = MKD_gg2aa.f):
C       RS      collider energy [GeV]                  (100000)
C       PDFSET  LHAPDF set name                        ('CT14lo')
C       MT, GT  top mass and width [GeV]               (173.0, 1.5)
C       MUG     scale of alpha_s in the Green function (40)
C       ETMAX   photon |eta| cut                       (2.5)
C       PTMIN   photon pT cut [GeV]                    (40)
C       MAAMIN, MAAMAX, DMAA  M_aa grid [GeV]          (300, 400, 0.1)
C       NCALL, ITMX  VEGAS calls / iterations          (50000, 6)
C       OUTFILE output file                            (Tab_<MT>_<GT>.dat)
C     Output: one line per M_aa,  "M_aa  dsigma/dM_aa[fb/GeV]  error",
C     in the format read by TMDP.read_template_from_file.
      PROGRAM MKTEMPLATE
      IMPLICIT NONE
      INTEGER I,NPT
      DOUBLE PRECISION MT,GT,MU,ASG,RS,MAA,ETMAX,PTMIN
      COMMON /SIGAA/   MT,GT,MU,ASG,RS,MAA,ETMAX,PTMIN
      DOUBLE PRECISION ALP,ASR,MUR,MUF
      COMMON /COUP/    ALP,ASR,MUR,MUF
      DOUBLE PRECISION S1,S2,S3,S4
      COMMON /RESULT/  S1,S2,S3,S4
      DOUBLE PRECISION INT2
      EXTERNAL         INT2
      DOUBLE PRECISION DSDMAA,ERR
      DOUBLE PRECISION MUG,MAAMIN,MAAMAX,DMAA
      INTEGER NCALL,ITMX
      CHARACTER*64  PDFSET
      CHARACTER*256 OUTFILE
      NAMELIST /TMPL/ RS,PDFSET,MT,GT,MUG,ETMAX,PTMIN,
     -     MAAMIN,MAAMAX,DMAA,NCALL,ITMX,OUTFILE
      INCLUDE 'parameter.inc'
      INCLUDE 'qcdparam.inc'
      INCLUDE 'qcdfunc.inc'
C.....Defaults (as in MKD_gg2aa.f)
      RS     = 100D3
      PDFSET = 'CT14lo'
      MT     = 173D0
      GT     = 1.5D0
      MUG    = 40D0
      ETMAX  = 2.5D0
      PTMIN  = 40D0
      MAAMIN = 300D0
      MAAMAX = 400D0
      DMAA   = 0.1D0
      NCALL  = 50000
      ITMX   = 6
      OUTFILE = ' '
      READ (5,NML=TMPL)
      IF ( OUTFILE.EQ.' ' )
     -     WRITE (OUTFILE,'("Tab_",F5.1,"_",F4.2,".dat")') MT,GT
C.....
      ALP = 1D0/128D0
      NQCD = 0
      CALL QCDLIBS
      MU  = MUG
      ASG = ASQCD(MU)
      CALL INITPDFSETBYNAME (PDFSET)
C.....
      WRITE (6,'(A,F9.1,A,A,A,F7.2,A,F5.2)') ' # RS=',RS,' PDF=',
     -     TRIM(PDFSET),' MT=',MT,' GT=',GT
      OPEN (10,FILE=OUTFILE,STATUS='REPLACE')
      NPT = NINT((MAAMAX-MAAMIN)/DMAA) + 1
      DO 100 I = 1, NPT
         MAA = MAAMIN + DMAA*(I-1)
         MUF = MAA
         MUR = MAA
         ASR = ASQCD (MUR)
         CALL VEGAS (INT2,1D-4,2,NCALL,ITMX,1,0)
         DSDMAA = S1 * 2D0*MAA/RS**2 ! [fb/GeV]
         ERR    = S2 * 2D0*MAA/RS**2
         WRITE (10,'(F12.4,1X,2(1PE15.5))') MAA, DSDMAA, ERR
 100  CONTINUE
      CLOSE (10)
      STOP
      END
C
C     Integrand: identical to INT2 in MKD_gg2aa.f.
      DOUBLE PRECISION FUNCTION INT2 (X)
      IMPLICIT NONE
      DOUBLE PRECISION X(2)
      DOUBLE PRECISION MT,GT,MU,ASG,RS,MAA,ETMAX,PTMIN
      COMMON /SIGAA/   MT,GT,MU,ASG,RS,MAA,ETMAX,PTMIN
      DOUBLE PRECISION ALP,ASR,MUR,MUF
      COMMON /COUP/    ALP,ASR,MUR,MUF
      DOUBLE PRECISION ET1,ET1MAX,ET1MIN,ET1JAC
      DOUBLE PRECISION ET2,ET2MAX,ET2MIN,ET2JAC
      DOUBLE PRECISION R,TAU,Y,ETHAT,COS,PTA
      DOUBLE PRECISION X1,X2
      DOUBLE COMPLEX   GG2AAQ,GG2AAT,GG2AAG
      EXTERNAL         GG2AAQ,GG2AAT,GG2AAG
      DOUBLE COMPLEX   AMP,AMPQ,AMPT,AMPG
      DOUBLE PRECISION SAMP1,SAMP2,SAMP3,SAMP4,SAMP5
      DOUBLE PRECISION SAMP,FAC,JAC
      DOUBLE PRECISION Q2Q,Q2T
      PARAMETER ( Q2Q = 1.22222222D0, Q2T = 0.4444444444D0 )
      DOUBLE PRECISION ONE,TWO,FOUR,HALF, PI
      PARAMETER ( ONE = 1D0, TWO = 2D0, FOUR = 4D0, HALF = 0.5D0 )
      PARAMETER ( PI = 3.141592654D0 )
      INTEGER I,L1234L(4,16)
      DATA (L1234L(I, 1),I=1,4) / 1, 1, 1, 1/
      DATA (L1234L(I, 2),I=1,4) / 1, 1, 1,-1/
      DATA (L1234L(I, 3),I=1,4) / 1, 1,-1, 1/
      DATA (L1234L(I, 4),I=1,4) / 1, 1,-1,-1/
      DATA (L1234L(I, 5),I=1,4) / 1,-1, 1, 1/
      DATA (L1234L(I, 6),I=1,4) / 1,-1, 1,-1/
      DATA (L1234L(I, 7),I=1,4) / 1,-1,-1, 1/
      DATA (L1234L(I, 8),I=1,4) / 1,-1,-1,-1/
      DATA (L1234L(I, 9),I=1,4) /-1, 1, 1, 1/
      DATA (L1234L(I,10),I=1,4) /-1, 1, 1,-1/
      DATA (L1234L(I,11),I=1,4) /-1, 1,-1, 1/
      DATA (L1234L(I,12),I=1,4) /-1, 1,-1,-1/
      DATA (L1234L(I,13),I=1,4) /-1,-1, 1, 1/
      DATA (L1234L(I,14),I=1,4) /-1,-1, 1,-1/
      DATA (L1234L(I,15),I=1,4) /-1,-1,-1, 1/
      DATA (L1234L(I,16),I=1,4) /-1,-1,-1,-1/
      DOUBLE PRECISION PDF1(-6:6), PDF2(-6:6)
      EXTERNAL EVOLVEPDF
C.....
      INT2 = 0D0
      R   = MAA/RS
      TAU = R**2
C...  Integrate over Eta1 and Eta2
      ET1MAX =  ETMAX
      ET1MIN = -ET1MAX
      ET1JAC = ET1MAX - ET1MIN
      ET1    = ET1JAC * X(1) + ET1MIN
      ET2MAX = MIN(-ET1 - DLOG(TAU), ETMAX)
      ET2MIN = MAX(-ET1 + DLOG(TAU),-ETMAX)
      ET2JAC = ET2MAX - ET2MIN
      ET2    = ET2JAC * X(2) + ET2MIN
      Y      = (ET1 + ET2) * HALF
      ETHAT  = (ET1 - ET2) * HALF
      COS    = DTANH(ETHAT)
      JAC    = HALF * ET1JAC * ET2JAC / DCOSH(ETHAT)**2
C...  Kinematical Cuts
      PTA = MAA*DSQRT(ONE-COS**2) * HALF
      IF ( DABS(ET1).GT.ETMAX ) RETURN
      IF ( DABS(ET2).GT.ETMAX ) RETURN
      IF ( PTA.LT.PTMIN ) RETURN
      IF ( PTA.LT.0.4*MAA ) RETURN
C...  Gluon Distribution Function
      X1  = R * DEXP( Y)
      X2  = R * DEXP(-Y)
      CALL EVOLVEPDF (X1,MUF,PDF1)
      CALL EVOLVEPDF (X2,MUF,PDF2)
C---- Matrix Elements Square ----
C...  PPPP
      AMPQ = GG2AAQ (   MAA,COS,L1234L(1,1))
      AMPT = GG2AAT (MT,MAA,COS,L1234L(1,1))
      AMPG = GG2AAG (MT,GT,MU,ASG,MAA,COS,L1234L(1,1),3)
      AMP  = Q2Q*AMPQ + Q2T*(AMPT+AMPG)
      SAMP1 = DBLE(AMP*DCONJG(AMP))
C...  MPPP
      AMPQ = GG2AAQ (   MAA,COS,L1234L(1,2))
      AMPT = GG2AAT (MT,MAA,COS,L1234L(1,2))
      AMPG = GG2AAG (MT,GT,MU,ASG,MAA,COS,L1234L(1,2),3)
      AMP  = Q2Q*AMPQ + Q2T*(AMPT+AMPG)
      SAMP2 = DBLE(AMP*DCONJG(AMP))
C...  MMPP
      AMPQ = GG2AAQ (   MAA,COS,L1234L(1,4))
      AMPT = GG2AAT (MT,MAA,COS,L1234L(1,4))
      AMPG = GG2AAG (MT,GT,MU,ASG,MAA,COS,L1234L(1,4),3)
      AMP  = Q2Q*AMPQ + Q2T*(AMPT+AMPG)
      SAMP3 = DBLE(AMP*DCONJG(AMP))
C...  MPMP
      AMPQ = GG2AAQ (   MAA,COS,L1234L(1,6))
      AMPT = GG2AAT (MT,MAA,COS,L1234L(1,6))
      AMPG = GG2AAG (MT,GT,MU,ASG,MAA,COS,L1234L(1,6),3)
      AMP  = Q2Q*AMPQ + Q2T*(AMPT+AMPG)
      SAMP4 = DBLE(AMP*DCONJG(AMP))
C...  MPPM
      AMPQ = GG2AAQ (   MAA,COS,L1234L(1,7))
      AMPT = GG2AAT (MT,MAA,COS,L1234L(1,7))
      AMPG = GG2AAG (MT,GT,MU,ASG,MAA,COS,L1234L(1,7),3)
      AMP  = Q2Q*AMPQ + Q2T*(AMPT+AMPG)
      SAMP5 = DBLE(AMP*DCONJG(AMP))
      SAMP = 2D0*SAMP1 + 8D0*SAMP2 + 2D0*SAMP3 + 2D0*SAMP4 + 2D0*SAMP5
C...  Hadronic Cross Section
      FAC  = ALP**2*ASR**2/(128D0*PI*MAA**2) * JAC * 389429.57D6 ! [fb]
      INT2 = FAC * SAMP * PDF1(0) / X1 * PDF2(0) / X2 ! dSigma / dTau [fb]
      RETURN
      END
