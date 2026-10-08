      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER I,J,IHEL,JAMP
      DOUBLE PRECISION MT,RS
      COMMON /PWD/     MT,RS,IHEL,JAMP
      DOUBLE COMPLEX INT1
      EXTERNAL       INT1
      DOUBLE PRECISION S1,S2,S3,S4
      COMMON /RESULT/  S1,S2,S3,S4
      DOUBLE PRECISION ANS(3),ERR(3)
      INTEGER L1234L(4,16)
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
      INCLUDE 'parameter.inc'
      RS = 300D0
      MT = RS/2D0
C.....
      IHEL = 2
C.....
      JAMP = 0
      CALL VEGAS (INT1,1D-5,1, 100000,6,0,0)
      ANS(1) = S1
      ERR(1) = S2
C.....
      JAMP = 2
      CALL VEGAS (INT1,1D-5,1, 500000,6,0,0)
      ANS(2) = S1
      ERR(2) = S2
C.....
      JAMP = 4
      CALL VEGAS (INT1,1D-5,1,2000000,6,0,0)
c     CALL BSINIT
c     CALL USERIN
c     CALL BASES (INT1,S1,S2,CTIME,IT1,IT2)
c     CALL BSINFO (6)
c     CALL BHPLOT (6)
      ANS(3) = S1
      ERR(3) = S2
C.....
      WRITE (6,'(4(I2,1X),1X,3(1PE15.4))') (L1234L(I,IHEL),I=1,4), ANS
      WRITE (6,'(4(I2,1X),1X,3(1PE15.4))') (L1234L(I,IHEL),I=1,4), ERR
C.....
      STOP
      END

      DOUBLE PRECISION FUNCTION INT1 (X)
      IMPLICIT NONE
      DOUBLE PRECISION X(1)
      DOUBLE PRECISION MT,RS
      INTEGER                IHEL,J
      COMMON /PWD/     MT,RS,IHEL,J
      INTEGER L1234L(4,16)
      COMMON /HEL/ L1234L
      DOUBLE PRECISION COS,CMIN,CMAX,CJAC
      INTEGER M1,M2
      DOUBLE COMPLEX   GG2AAT
      DOUBLE PRECISION AMP,AMP1,GG2AAA,WIGNERD
      EXTERNAL         GG2AAT,  GG2AAA
      CMIN = -1D0
      CMAX =  1D0
      CJAC = CMAX - CMIN
      COS  = CJAC*X(1) + CMIN
      M1   = ABS(L1234L(2,IHEL) - L1234L(1,IHEL))
      M2   = ABS(L1234L(3,IHEL) - L1234L(4,IHEL))
c     AMP  = DREAL(GG2AAT (MT,RS,COS,L1234L(1,IHEL)))
      AMP  = GG2AAA(COS,L1234L(1,IHEL))
      INT1 = AMP * WIGNERD (J,M1,M2,COS) / 2D0 * CJAC
      RETURN
      END
      
      DOUBLE PRECISION FUNCTION WIGNERD (J,M1,M2,COS)
      IMPLICIT NONE
      INTEGER J,M1,M2
      DOUBLE PRECISION COS,DFN(0:4,0:4,0:4)
      DFN(0,0,0) =                          1D0
      DFN(2,0,0) = (             3D0*COS**2-1D0) / 2D0
      DFN(4,0,0) = (35D0*COS**4-30D0*COS**2+3D0) / 8D0
      DFN(2,2,0) = (1D0-COS**2) * DSQRT(3D0/8D0)
      DFN(4,2,0) = (1D0-COS**2) * (7D0*COS**2-1D0) * DSQRT(5D0/32D0)
      DFN(2,2,2) = (1D0+COS)**2 / 4D0
      DFN(3,2,2) = (1D0+COS)**2 / 4D0 * (3D0*COS-2D0)
      DFN(4,2,2) = (1D0+COS)**2 / 4D0 * (7D0*COS**2-7D0*COS+1D0)
      DFN(2,0,2) = DFN(2,2,0)
      DFN(4,0,2) = DFN(4,2,0)
      WIGNERD    = DFN(J,M1,M2)
      RETURN
      END
      
