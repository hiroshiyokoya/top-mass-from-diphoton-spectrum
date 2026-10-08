      DOUBLE COMPLEX FUNCTION GG2AAT (MT,RS,COS,L1234)
      IMPLICIT NONE
      DOUBLE PRECISION MT,RS,COS
      INTEGER L1234(4)
      DOUBLE PRECISION ONE,TWO,FOUR,HALF
      PARAMETER ( ONE = 1D0, TWO = 2D0, FOUR = 4D0, HALF = 0.5D0 )
      DOUBLE PRECISION S,T,U,R
      DOUBLE COMPLEX MTPPPP,MTMPPP,MTMMPP
      DOUBLE COMPLEX BFN,TFN,IFN
      EXTERNAL       BFN,TFN,IFN
      MTPPPP (R,S,T) = ONE - HALF/(R*S)*IFN(R,S)
     -     - HALF/(R*T)*IFN(R,T) - HALF/(S*T)*IFN(S,T)
      MTMPPP (R,S,T) = ONE + (ONE/R+ONE/S+ONE/T)
     -     * (TFN(R)+TFN(S)+TFN(T))
     -     - (ONE/T+HALF/(R*S))*IFN(R,S)
     -     - (ONE/S+HALF/(R*T))*IFN(R,T)
     -     - (ONE/R+HALF/(S*T))*IFN(S,T)
      MTMMPP (R,S,T) = - ONE - (TWO+FOUR*S/R)*BFN(S)
     -     - (TWO+FOUR*T/R)*BFN(T)
     -     - TWO*((S**2+T**2)/R**2-ONE/R) * (TFN(S)+TFN(T))
     -     + (ONE/S-HALF/(R*S))*IFN(R,S)
     -     + (ONE/T-HALF/(R*T))*IFN(R,T)
     -     + (TWO*(S**2+T**2)/R**2-FOUR/R-ONE/S-ONE/T-HALF/(S*T))*IFN(S,T)
c     IF ( DABS(COS-ONE) .LE. 1D-4 ) COS =  1D0 - 1D-2
c     IF ( DABS(COS+ONE) .LE. 1D-4 ) COS = -1D0 + 1D-2
      S = RS**2
      T = - HALF*S*(ONE-COS)
      U = - HALF*S*(ONE+COS)
      R = S/(4D0*MT**2)
      S = T/(4D0*MT**2)
      T = U/(4D0*MT**2)
      IF     ( ABS(L1234(1)+L1234(2)+L1234(3)+L1234(4)).EQ.4 ) THEN
         GG2AAT = MTPPPP (R,S,T)
      ELSEIF ( ABS(L1234(1)+L1234(2)+L1234(3)+L1234(4)).EQ.2 ) THEN
         GG2AAT = MTMPPP (R,S,T)
      ELSEIF ( ABS(L1234(1)+L1234(2)-L1234(3)-L1234(4)).EQ.4 ) THEN
         GG2AAT = MTMMPP (R,S,T)
      ELSEIF ( ABS(L1234(1)-L1234(2)-L1234(3)+L1234(4)).EQ.4 ) THEN
         GG2AAT = MTMMPP (T,S,R)
      ELSEIF ( ABS(L1234(1)-L1234(2)+L1234(3)-L1234(4)).EQ.4 ) THEN
         GG2AAT = MTMMPP (S,T,R)
      ELSE
         STOP
      ENDIF
      RETURN
      END
      
      DOUBLE COMPLEX FUNCTION BFN (R)
      IMPLICIT NONE
      DOUBLE PRECISION R,B
      DOUBLE COMPLEX L
      DOUBLE PRECISION ZERO,ONE,HALF,PI
      PARAMETER ( ZERO = 0D0, ONE = 1D0, HALF = 0.5D0 )
      PARAMETER ( PI = 3.141592654D0 )
      DOUBLE COMPLEX IMAG
      PARAMETER ( IMAG = DCMPLX(0D0,1D0) )
      IF ( R.GE.ONE ) THEN
         B = DSQRT(ONE-ONE/R)
         BFN = B * DACOSH(DSQRT(R)) - ONE - IMAG*PI*HALF*B
      ELSE IF ( R.GT.ZERO ) THEN
         B = DSQRT(ONE/R-ONE)
         BFN = B * DASIN(DSQRT(R)) - ONE
      ELSE IF ( R.EQ.ZERO ) THEN
         BFN = ZERO
      ELSE IF ( R.LT.ZERO ) THEN
         B = DSQRT(ONE-ONE/R)
         BFN = B * DASINH(DSQRT(-R)) - ONE
      ENDIF
      RETURN
      END
      
      DOUBLE COMPLEX FUNCTION TFN (R)
      IMPLICIT NONE
      DOUBLE PRECISION R,B
      DOUBLE COMPLEX L
      DOUBLE PRECISION ZERO,ONE,HALF,PI
      PARAMETER ( ZERO = 0D0, ONE = 1D0, HALF = 0.5D0 )
      PARAMETER ( PI = 3.141592654D0 )
      DOUBLE COMPLEX IMAG
      PARAMETER ( IMAG = DCMPLX(0D0,1D0) )
      IF ( R.GE.ONE ) THEN
         TFN = DACOSH(DSQRT(R))**2 - PI**2/4D0 - IMAG*PI*DACOSH(DSQRT(R))
      ELSE IF ( R.GE.ZERO ) THEN
         TFN = - DASIN(DSQRT(R))**2
      ELSE IF ( R.LT.ZERO ) THEN
         TFN = DASINH(DSQRT(-R))**2
      ENDIF
      RETURN
      END
      
      DOUBLE COMPLEX FUNCTION IFN (R,S)
      IMPLICIT NONE
      DOUBLE PRECISION R,S,A
      DOUBLE COMPLEX B,DLA,DLB,DLC,DLD
      DOUBLE PRECISION ZERO,ONE,TWO,PI
      PARAMETER ( ZERO = 0D0, ONE = 1D0, TWO = 2D0 )
      PARAMETER ( PI = 3.141592654D0 )
      DOUBLE COMPLEX IMAG
      PARAMETER ( IMAG = DCMPLX(0D0,1D0) )
      INCLUDE 'CHAPLIN.INC'
      double precision ddilog
      A   = DSQRT( ONE - (R+S)/(R*S) )
      IF ( R.LT.ZERO .OR. R.GT.ONE ) THEN
         B = DSQRT(ONE-ONE/R)
      ELSE
         B = IMAG*DSQRT(ONE/R-ONE)
      ENDIF
      Z = (A+ONE)/(A+B)
      DLA = HPL2(0,1,Z)
      Z = (A-ONE)/(A+B)
      DLB = HPL2(0,1,Z)
      Z = (A+ONE)/(A-B)
      DLC = HPL2(0,1,Z)
      Z = (A-ONE)/(A-B)
      DLD = HPL2(0,1,Z)
      IFN = DREAL( - DLA + DLB - DLC + DLD )/(2D0*A)
      IF ( R.GT.ONE ) IFN = IFN + IMAG*PI/(TWO*A)*LOG((A-B)/(A+B))
      IF ( S.LT.ZERO .OR. S.GT.ONE ) THEN
         B = DSQRT(ONE-ONE/S)
      ELSE
         B = IMAG*DSQRT(ONE/S-ONE)
      ENDIF
      Z = (A+ONE)/(A+B)
      DLA = HPL2(0,1,Z)
      Z = (A-ONE)/(A+B)
      DLB = HPL2(0,1,Z)
      Z = (A+ONE)/(A-B)
      DLC = HPL2(0,1,Z)
      Z = (A-ONE)/(A-B)
      DLD = HPL2(0,1,Z)
      IFN = IFN + DREAL( - DLA + DLB - DLC + DLD )/(2D0*A)
      IF ( S.GT.ONE ) IFN = IFN + IMAG*PI/(TWO*A)*LOG((A-B)/(A+B))
      RETURN
      END
