      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER I,J,IHEL,JAMP
      DOUBLE PRECISION RS,MAA,ETMAX,PTMIN
      COMMON /SIGAA/   RS,MAA,ETMAX,PTMIN
      DOUBLE PRECISION ALP,ASR,MUR
      COMMON /COUP/    ALP,ASR,MUR
      DOUBLE COMPLEX INT2
      EXTERNAL       INT2
      DOUBLE PRECISION S1,S2,S3,S4
      COMMON /RESULT/  S1,S2,S3,S4
      DOUBLE PRECISION DSDMAA, ERR
      INCLUDE 'parameter.inc'
      INCLUDE 'qcdparam.inc'
      INCLUDE 'qcdfunc.inc'
      DOUBLE PRECISION MTT(105),SIG(105),DER,H
      INTEGER NB
C.....
      ALP = 1D0/128D0
c     NQCD = 0
c     CALL QCDLIBS
C.....
      RS = 13D3
      ETMAX = 2.5D0
      PTMIN =  40D0
C...  LHAPDF Initialization
      CALL INITPDFSETBYNAME ('CT14nlo')
C.....
      H = 5D0
      NB = 100/INT(H)
      DO 100 I = 1, NB+5
         MAA = 300D0 + H*(I-3)
c     MAA = 350D0
C.....
c     MUR = MAA
c     ASR = ASQCD (MUR)
c     PRINT *, ALP,ASR,MUR
C.....
      CALL VEGAS (INT2,1D-4,2,10000,6,0,0)
C.....
c         CALL BSINIT
c         CALL USERIN
c         CALL BASES (INT2,S1,S2,CTIME,IT1,IT2)
c         CALL BSINFO (6)
c         CALL BSINFO (7)
c         CALL BHPLOT (7)
C.....
         DSDMAA = S1 * 2D0*MAA/RS**2 ! [fb/GeV]
         ERR    = S2 * 2D0*MAA/RS**2
         WRITE ( 6,'(F12.4,1X,2(1PE15.5)))') MAA, DSDMAA, ERR
         WRITE (11,'(F12.4,1X,2(1PE15.5)))') MAA, DSDMAA, ERR
C.....
         MTT(I) = MAA
         SIG(I) = DSDMAA
C.....
 100  CONTINUE
C.....
      DO 101 I = 3, NB+3
         MAA = MTT(I)
         DER = (SIG(I-2)-8D0*SIG(I-1)+8D0*SIG(I+1)-SIG(I+2)) / (12D0*H)
         WRITE ( 6,'(F12.4,1X,1(1PE15.5)))') MAA, DER/SIG(I)*MAA
         WRITE (21,'(F12.4,1X,1(1PE15.5)))') MAA, DER/SIG(I)*MAA
 101  CONTINUE
C.....
      STOP
      END
      
      DOUBLE PRECISION FUNCTION INT2 (X)
      IMPLICIT NONE
      DOUBLE PRECISION X(2)
      DOUBLE PRECISION RS,MAA,ETMAX,PTMIN
      COMMON /SIGAA/   RS,MAA,ETMAX,PTMIN
      DOUBLE PRECISION ALP,ASR,MUR
      COMMON /COUP/    ALP,ASR,MUR
      DOUBLE PRECISION ET1,ET1MAX,ET1MIN,ET1JAC
      DOUBLE PRECISION ET2,ET2MAX,ET2MIN,ET2JAC
      DOUBLE PRECISION YMAX,YMIN,YJAC
      DOUBLE PRECISION CMAX,CMIN,CJAC
      DOUBLE PRECISION R,TAU,Y,ETHAT,COS,PTA
      DOUBLE PRECISION X1,X2,MUF,LUM
      DOUBLE PRECISION SAMP,FAC,JAC
      DOUBLE PRECISION ONE,TWO,FOUR,HALF, PI
      PARAMETER ( ONE = 1D0, TWO = 2D0, FOUR = 4D0, HALF = 0.5D0 )
      PARAMETER ( PI = 3.141592654D0 )
      INTEGER I
C.....
      DOUBLE PRECISION PDF1(-6:6), PDF2(-6:6)
      EXTERNAL EVOLVEPDF
C.....
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
C...  Integrate over Y and Cos
c     YMAX = -DLOG(R)
c     YMIN =  DLOG(R)
c     YJAC = YMAX - YMIN
c     Y    = YJAC*X(1) + YMIN
c     CMAX = 1D0
c     CMIN = -CMAX
c     CJAC = CMAX - CMIN
c     COS  = CJAC*X(2) + CMIN
c     JAC  = YJAC * CJAC
c     ETHAT = HALF * DLOG((ONE+COS)/(ONE-COS))
c     ET1 =  ETHAT + Y
c     ET2 = -ETHAT + Y
C...  Kinematical Cuts
      PTA = MAA*DSQRT(ONE-COS**2) * HALF
      IF ( DABS(ET1).GT.ETMAX ) RETURN
      IF ( DABS(ET2).GT.ETMAX ) RETURN
      IF ( PTA.LT.PTMIN ) RETURN
C...  Gluon Distribution Function
      MUF = MAA
      X1  = R * DEXP( Y)
      X2  = R * DEXP(-Y)
      CALL EVOLVEPDF (X1,MUF,PDF1)
      CALL EVOLVEPDF (X2,MUF,PDF2)
      LUM =  (4D0/9D0)**2 * ( PDF1( 2)*PDF2(-2) + PDF1(-2)*PDF2( 2)
     -     +                  PDF1( 4)*PDF2(-4) + PDF1(-4)*PDF2( 4) )
     -     + (1D0/9D0)**2 * ( PDF1( 1)*PDF2(-1) + PDF1(-1)*PDF2( 1)
     -     +                  PDF1( 3)*PDF2(-3) + PDF1(-3)*PDF2( 3)
     -     +                  PDF1( 5)*PDF2(-5) + PDF1(-5)*PDF2( 5) )
C---- Matrix Elements Square
      SAMP = DCOSH(2D0*ETHAT)
C...  Hadronic Cross Section
      FAC  = PI*ALP**2/(3D0*MAA**2) * JAC * 389429.57D6 ! [fb]
      INT2 = FAC * SAMP * LUM / (X1*X2) ! dSigma / dTau [fb]
C.....
      CALL XHFILL (1,X(1),INT2)
      CALL XHFILL (2,X(2),INT2)
      CALL XHFILL (3,ET1 ,INT2)
      CALL XHFILL (4,ET2 ,INT2)
C.....
      RETURN
      END

      SUBROUTINE USERIN
      IMPLICIT NONE
      INCLUDE 'parameter.inc'
      INTEGER MAXDIM
      PARAMETER (MAXDIM=100)
      DOUBLE PRECISION XL(MAXDIM),XU(MAXDIM)
      INTEGER IG(MAXDIM)
      INTEGER NCALL, ITMX1, ITMX2
      DOUBLE PRECISION ACC1, ACC2
      INTEGER IDIM
C.....
      DOUBLE PRECISION RS,MAA,ETMAX,PTMIN
      COMMON /SIGAA/   RS,MAA,ETMAX,PTMIN
C
      NDIM  = 2
      NWILD = 2
C
      DO IDIM = 1, NDIM
         XL(IDIM) = 0D0
         XU(IDIM) = 1D0
         IG(IDIM) = 1
      ENDDO
C
      CALL BSSETD (NDIM,NWILD,XL,XU,IG)
C
      NCALL = 10 000
      ITMX1 = 6
      ITMX2 = 6
      ACC1  = 0.1D0
      ACC2  = 0.1D0
C
      CALL BSSETP (NCALL,ITMX1,ITMX2,ACC1,ACC2)
      CALL XHINIT ( 1, 0D0,1D0,50,'D SIGMA / D X(1)')
      CALL XHINIT ( 2, 0D0,1D0,50,'D SIGMA / D X(2)')
      CALL XHINIT ( 3,-ETMAX,ETMAX,50,'D SIGMA / D Et1')
      CALL XHINIT ( 4,-ETMAX,ETMAX,50,'D SIGMA / D Et2')
C
      RETURN
      END
