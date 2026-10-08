      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER I,J,IHEL,JAMP,ICMP
      DOUBLE PRECISION MT,GT,MU,ASG,RS
      COMMON /PWD/     MT,GT,MU,ASG,RS,IHEL,JAMP,ICMP
      DOUBLE COMPLEX INT1
      EXTERNAL       INT1
      DOUBLE PRECISION S1,S2,S3,S4
      COMMON /RESULT/  S1,S2,S3,S4
      DOUBLE PRECISION ReM,ImM
      DOUBLE PRECISION Req,Imq,Q2Q,Q2T
      DOUBLE COMPLEX GG2AAG,GG2AAT
      EXTERNAL       GG2AAG,GG2AAT
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
      INCLUDE 'qcdparam.inc'
      INCLUDE 'qcdfunc.inc'
      CALL QCDLIBS
      MT = 173D0
      GT = 1.498D0
      IHEL = 13
      JAMP = 0
C.....
      REQ = -1./3. - 4.*PI**2/9.
      IMQ = 0.
      Q2Q = 11./9.
      Q2T = 4./9.
C.....
      MU  = 40D0
      ASG = ASQCD(MU)
C.....
      DO I = 1,41
         RS = 342D0 + 0.1D0*(I-1)
c      DO I = 1,301
c         RS = 330D0 + 0.1D0*(I-1)
c     RS = 360D0 + .01D0*1.1**(I-1)
C.....
         ICMP  = 0
         CALL VEGAS (INT1,1D-3,1, 10000,6,0,0)
         ReM = S1
C.....
c         IF ( RS.GT.2D0*MT ) THEN
c            ICMP  = 1
c            CALL VEGAS (INT1,1D-3,1, 10000,6,0,0)
c            ImM = S1
c         ELSE
c            ImM = 0D0
c         ENDIF
c         ImM = DIMAG( GG2AAT(MT,RS,0D0,L1234L(1,IHEL)) )
c         ImM = DIMAG( GG2AAT(MT,RS,0D0,L1234L(1,IHEL)) + GG2AAG(MT,GT,MU,ASG,RS,0D0,L1234L(1,IHEL),2) )
         ImM = DIMAG( GG2AAG(MT,GT,MU,ASG,RS,0D0,L1234L(1,IHEL),2) )
c         ENDIF
C.....
         WRITE ( 6,*) RS, -ReM,-ImM
         WRITE (11,*) RS, -Q2T*ReM, -Q2T*ImM
         WRITE (12,*) RS, -Q2Q*Req - Q2T*ReM, -Q2Q*Imq - Q2T*ImM
C.....
      ENDDO
C.....
      STOP
      END

      DOUBLE PRECISION FUNCTION INT1 (X)
      IMPLICIT NONE
      DOUBLE PRECISION X(1)
      DOUBLE PRECISION MT,GT,MU,ASG,RS
      INTEGER                IHEL,J,ICMP
      COMMON /PWD/     MT,GT,MU,ASG,RS,IHEL,J,ICMP
      INTEGER L1234L(4,16)
      COMMON /HEL/ L1234L
      DOUBLE PRECISION COS,CMIN,CMAX,CJAC
      INTEGER M1,M2
      DOUBLE COMPLEX   CAMP,GG2AAT,GG2AAG
      DOUBLE PRECISION  AMP,GG2AAA,WIGNERD
      EXTERNAL              GG2AAT,GG2AAG,GG2AAA
      CMIN = -1D0
      CMAX =  1D0
      CJAC = CMAX - CMIN
      COS  = CJAC*X(1) + CMIN
      M1   = ABS(L1234L(2,IHEL) - L1234L(1,IHEL))
      M2   = ABS(L1234L(3,IHEL) - L1234L(4,IHEL))
      CAMP = GG2AAG (MT,GT,MU,ASG,RS,COS,L1234L(1,IHEL),2)
      IF ( ICMP.EQ.0 ) THEN
         AMP = DREAL(CAMP)
      ELSE
         AMP = DIMAG(CAMP)
      ENDIF
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
