      SUBROUTINE GRNLOMSB (E,M,GAM,I,AS,MU,REG,IMG)
C---  I 0:G0 with no width, 1:O(as) with no width, 2:G0 with finite width,
C---  I 3:O(as) with finite width, 4:LO
      IMPLICIT NONE
      INTEGER I,N
      DOUBLE PRECISION E,M,GAM,AS,MU,REG,IMG
      DOUBLE PRECISION G,P0,P1,P2,PN,D
      DOUBLE PRECISION PI
      PARAMETER ( PI=3.14159265359D0 )
      Double precision Y0S,EB
      double complex BW
      REG = 0D0
      IMG = 0D0
      IF ( I.LE.1 ) THEN
         G = 1D-12 * M
      ELSE
         G = GAM
      ENDIF
      P1 = DSQRT(0.5D0*M*(DSQRT(E**2+G**2)+E))
      P2 = DSQRT(0.5D0*M*(DSQRT(E**2+G**2)-E))
C.....
      REG = -M*P2/(4D0*PI)
      IMG =  M*P1/(4D0*PI)
      IF ( I.EQ.0 .OR. I.EQ.2 ) RETURN
C.....
      P0 = 2D0/3D0*M*AS
      D  = MU/(2D0*M) * DEXP(0.5D0) ! MSB
      REG = REG + M*P0/(4D0*PI) * DLOG(M**2/(P1**2+P2**2)*D**2)
      IMG = IMG + M*P0/(2D0*PI) * DATAN2(P1,P2)
      IF ( I.LE.3 ) RETURN
C.....
      DO N = 1,1!300
         PN  = P0/N
c         REG = REG + M*P0**2/(2D0*PI)*(P2-PN)/N**2/((P2-PN)**2+P1**2)
c         IMG = IMG + M*P0**2/(2D0*PI)*P1/N**2/((P2-PN)**2+P1**2)
c         REG =  M*P0**2/(2D0*PI)*(P2-PN)/N**2/((P2-PN)**2+P1**2)
c         IMG =  M*P0**2/(2D0*PI)*P1/N**2/((P2-PN)**2+P1**2)
         Y0S = P0**3/PI/N**3
         EB  = -P0**2/M/N**2
         BW = - Y0S/DCMPLX(E-EB,GAM)
         REG = DREAL(BW)
         IMG = DIMAG(BW)
      ENDDO
      RETURN
      END
