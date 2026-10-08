      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER I,J,IHEL,JAMP
      DOUBLE PRECISION MT,MAA
      COMMON /SIGAA/   MT,MAA
      DOUBLE COMPLEX INT1
      EXTERNAL       INT1
      DOUBLE PRECISION S1,S2,S3,S4
      COMMON /RESULT/  S1,S2,S3,S4
      INCLUDE 'parameter.inc'
C.....
      MT = 173D0
C.....
      MAA = 350D0
C.....
      CALL VEGAS (INT1,1D-3,1,100000,6,0,0)
C.....
c     CALL BSINIT
c     CALL USERIN
c     CALL BASES (INT1,S1,S2,CTIME,IT1,IT2)
c     CALL BSINFO (6)
c     CALL BSINFO (7)
c     CALL BHPLOT (7)
C.....
      WRITE (6,'(F12.4,1X,2(1PE15.5)))') MAA, S1,S2
C.....
 100  CONTINUE
C.....
      STOP
      END

      DOUBLE PRECISION FUNCTION INT1 (X)
      IMPLICIT NONE
      DOUBLE PRECISION X(1)
      DOUBLE PRECISION MT,MAA
      COMMON /SIGAA/   MT,MAA
      DOUBLE PRECISION COS,CMAX,CMIN,CJAC
      DOUBLE PRECISION ET,ETMAX,ETMIN,ETJAC
      DOUBLE COMPLEX   GG2AAQ,GG2AAT
      EXTERNAL         GG2AAQ,GG2AAT
      DOUBLE COMPLEX   AMPQ,AMPT,AMP
      DOUBLE PRECISION JAC,SAMP
      DOUBLE PRECISION ONE,TWO,FOUR,HALF, PI
      PARAMETER ( ONE = 1D0, TWO = 2D0, FOUR = 4D0, HALF = 0.5D0 )
      PARAMETER ( PI = 3.141592654D0 )
      INTEGER I,IHEL,L1234L(4,16)
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
      COMMON /HEL/ L1234L
C.....
      CMAX =  1D0
      CMIN = -1D0
      CJAC = CMAX - CMIN
      COS  = CJAC*X(1) + CMIN
      JAC  = CJAC
C.....
c      ETMAX =  15D0
c      ETMIN = -15D0
c      ETJAC = ETMAX - ETMIN
c      ET    = ETJAC*X(1) + ETMIN
c      JAC   = ETJAC / DCOSH(ET)**2
c      COS   = DTANH(ET)
C...  Matrix Elements Squared
      IHEL = 1
      AMPQ = GG2AAQ (   MAA,COS,L1234L(1,IHEL))
c     AMPT = GG2AAT (MT,MAA,COS,L1234L(1,IHEL))
      SAMP = DBLE(AMPQ*DCONJG(AMPQ))
C...  Hadronic Cross Section
      INT1 = SAMP * JAC
C.....
      RETURN
      END
