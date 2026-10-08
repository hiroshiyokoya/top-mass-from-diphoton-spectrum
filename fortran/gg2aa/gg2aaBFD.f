      DOUBLE COMPLEX FUNCTION GG2AABFD1 (RS,COS,L1234)
      IMPLICIT NONE
      DOUBLE PRECISION RS,COS
      INTEGER L1234(4)
      DOUBLE PRECISION PI,ONE,HALF,EPS
      DOUBLE COMPLEX IMAG
      PARAMETER ( EPS = 1D-5 , PI = 3.141592654D0 )
      PARAMETER ( ONE = 1D0, HALF = 0.5D0, IMAG = DCMPLX(0D0,1D0) )
      DOUBLE PRECISION S,T,U
      DOUBLE COMPLEX MQPPPP,MQMPPP,MQMMPP,MQMPMP
      DATA MQPPPP /1D0/, MQMPPP /1D0/
      MQMMPP (S,T,U) = - HALF*(T**2+U**2)/S**2 * ( DLOG(T/U)**2
     -     + PI**2 ) - (T-U)/S * DLOG(T/U) - ONE
      MQMPMP (S,T,U) = - HALF*(T**2+S**2)/U**2 * LOG(-T/S)**2
     -     - (T-S)/U * LOG(-T/S) - ONE - IMAG*PI * (
     -     (T**2+S**2)/U**2 * LOG(-T/S) + (T-S)/U )
c     IF ( COS.EQ. 1D0 ) COS =  1D0-EPS
c     IF ( COS.EQ.-1D0 ) COS = -1D0+EPS
      S = RS**2
      T = - HALF*S*(ONE+COS)
      U = - HALF*S*(ONE-COS)
      IF     ( ABS(L1234(1)+L1234(2)+L1234(3)+L1234(4)).EQ.4 ) THEN
         GG2AABFD1 = MQPPPP
      ELSEIF ( ABS(L1234(1)+L1234(2)+L1234(3)+L1234(4)).EQ.2 ) THEN
         GG2AABFD1 = MQMPPP
      ELSEIF ( ABS(L1234(1)+L1234(2)-L1234(3)-L1234(4)).EQ.4 ) THEN
         GG2AABFD1 = MQMMPP (S,T,U)
      ELSEIF ( ABS(L1234(1)-L1234(2)+L1234(3)-L1234(4)).EQ.4 ) THEN
         GG2AABFD1 = MQMPMP (S,T,U)
      ELSEIF ( ABS(L1234(1)-L1234(2)-L1234(3)+L1234(4)).EQ.4 ) THEN
         GG2AABFD1 = MQMPMP (S,U,T)
      ELSE
         PRINT *, "Error: GG2AABFD1", L1234
         STOP
      ENDIF
      RETURN
      END
C
      DOUBLE COMPLEX FUNCTION GG2AABFD2 (RS,COS,MUR,L1234)
      IMPLICIT NONE
      DOUBLE PRECISION RS,COS,MUR
      INTEGER L1234(4)
      DOUBLE PRECISION PI,ONE,HALF,EPS
      DOUBLE COMPLEX IMAG
      INTEGER NC, NF
      PARAMETER ( NC = 3, NF = 5 )
      PARAMETER ( EPS = 1D-5 , PI = 3.141592654D0 )
      PARAMETER ( ONE = 1D0, HALF = 0.5D0, IMAG = DCMPLX(0D0,1D0) )
      DOUBLE PRECISION S,T,U
      DOUBLE COMPLEX AMP1, FL, FSL
      DOUBLE COMPLEX GG2AABFD1, BFDFL, BFDFSL
      EXTERNAL GG2AABFD1, BFDFL, BFDFSL
c     IF ( COS.EQ. 1D0 ) COS =  1D0-EPS
c     IF ( COS.EQ.-1D0 ) COS = -1D0+EPS
      S = RS**2
      T = - HALF*S*(ONE+COS)
      U = - HALF*S*(ONE-COS)
      AMP1 = GG2AABFD1(RS,COS,L1234)
      FL = BFDFL(S,T,U,L1234)
      FSL = BFDFSL(S,T,U,L1234)
      GG2AABFD2 = (11D0*NC-2D0*NF)/6D0 * (DLOG(MUR**2/S)+IMAG*PI) * AMP1
     -     + NC*FL - 1D0/NC*FSL
      RETURN
      END
C
      DOUBLE COMPLEX FUNCTION BFDFL (S,T,U,L1234)
      IMPLICIT NONE
      DOUBLE PRECISION S,T,U
      INTEGER L1234(4)
      DOUBLE PRECISION PI,ZETA3,ZETA4,ONE,HALF,TWO,EPS
      DOUBLE COMPLEX UNIT,IMAG
      PARAMETER ( EPS = 1D-5 , PI = 3.141592654D0 )
      PARAMETER ( ZETA3 = 1.2020569D0 , ZETA4 = PI**2/90D0 )
      PARAMETER ( ONE = 1D0, HALF = 0.5D0, TWO = 2D0 )
      PARAMETER ( UNIT = DCMPLX(1D0,0D0), IMAG = DCMPLX(0D0,1D0) )
      DOUBLE PRECISION X,Y,XX,YY
      DOUBLE COMPLEX POLYLOG
      EXTERNAL POLYLOG
      DOUBLE COMPLEX MPPP,PPMP,MMPP,MPMP
      MPPP(X,Y,XX,YY) = UNIT/8D0 * (
     -     (2D0+4D0*X/Y**2-5D0*X**2/Y**2) * ((XX+IMAG*PI)**2+PI**2)
     -     - (ONE-X*Y)*((XX-YY)**2+PI**2)
     -     + TWO*(9D0/Y-10D0*X)*(X+IMAG*PI) )
      PPMP(X,Y,XX,YY) = UNIT/8D0 * (
     -     (2D0+6D0*X/Y**2-3D0*X**2/Y**2) * ((XX+IMAG*PI)**2+PI**2)
     -     - (X-Y)**2*((XX-YY)**2+PI**2)
     -     + TWO*(9d0/Y- 8D0*X)*(X+IMAG*PI) )
      MMPP(X,Y,XX,YY) = -(X**2+Y**2) * (
     -     4D0*POLYLOG(4,-X) + (YY-3D0*XX-TWO*IMAG*PI)*POLYLOG(3,-X)
     -     + ((XX+IMAG*PI)**2+PI**2)*POLYLOG(2,-X) + (XX+YY)**4/48D0
     -     + IMAG*PI/12D0*(XX+YY)**3 + IMAG*HALF*PI**3*XX
     -     - PI**2/12D0*XX**2 - 109D0/720D0*PI**4 )
     -     + HALF*X*(ONE-3D0*Y) * (
     -     POLYLOG(3,-X/Y) - (XX-YY)*POLYLOG(2,-X/Y) - ZETA3
     -     + HALF*YY*((XX-YY)**2+PI**2) )
     -     + X**2/4D0 * (
     -     (XX-YY)**3 + 3D0*(YY+IMAG*PI)*((XX-YY)**2+PI**2) )
     -     + ONE/8D0*(14D0*(X-Y)-8D0/Y+9D0/Y**2)*((XX+IMAG*PI)**2+PI**2)
     -     + ONE/16D0*(38D0*X*Y-13D0)*((XX-YY)**2+PI**2)
     -     - PI**2/12D0 - 9d0/4D0*(ONE/Y+TWO*X)*(XX+IMAG*PI)
     -     + ONE/4D0
      MPMP(X,Y,XX,YY) = -TWO*(X**2+ONE)/Y**2 * ( POLYLOG(4,-X)
     -     - ZETA4 - HALF*(XX+IMAG*PI)*(POLYLOG(3,-X)-ZETA3)
     -     + PI**2/6D0*(POLYLOG(2,-X)-PI**2/6D0-HALF*XX**2) - XX**4/48D0
     -     + ONE/24D0*(XX+IMAG*PI)**2*((XX+IMAG*PI)**2+PI**2) )
     -     + TWO*(3D0*(ONE-X)**2-TWO)/Y**2 * (
     -     POLYLOG(4,-X) + POLYLOG(4,-X/Y) - POLYLOG(4,-Y)
     -     - (YY+IMAG*PI)*(POLYLOG(3,-X)-ZETA3)
     -     + PI**2/6D0*(POLYLOG(2,-X)+HALF*YY**2)
     -     - XX*YY**3/6D0 + YY**4/24D0 - 7D0/360D0*PI**4 )
     -     - TWO/3D0*(8D0-X+30D0*X/Y) * (POLYLOG(3,-Y) - ZETA3
     -     - (YY+IMAG*PI)*(POLYLOG(2,-Y)-PI**2/6D0)
     -     - HALF*XX*((YY+IMAG*PI)**2+PI**2))
     -     + ONE/6D0*(4D0*Y+27D0+42D0/Y+4D0/Y**2) * (
     -     POLYLOG(3,-X) - ZETA3 - (XX+IMAG*PI)*(POLYLOG(2,-X)-PI**2/6)
     -     + IMAG*HALF*PI*XX**2 - PI**2*XX )
     -     + ONE/12D0*(3D0-TWO/Y-12D0*X/Y**2)*(XX+IMAG*PI)
     -     * ((XX+IMAG*PI)**2+PI**2)
     -     - Y/3D0*(XX+IMAG*PI)*((YY+IMAG*PI)**2+PI**2)
     -     + TWO*(ONE+TWO/Y)*(ZETA3-PI**2/6D0*(YY+IMAG*PI))
     -     + 1./24D0*(Y**2-24D0*Y+44D0-8D0*X**3/Y)*((XX-YY)**2+PI**2)
     -     - 1./24D0*(15D0-14D0*X/Y-48D0*X/Y**2)*((XX+IMAG*PI)**2+PI**2)
     -     + 1./24D0*(8D0*X/Y+60D0-24D0*Y/X+27D0*Y**2/X**2)
     -     * ((YY+IMAG*PI)**2+PI**2)
     -     + 4D0/9D0*PI**2*X/Y
     -     + ONE/12D0*(2D0*X**2-54D0*X-27D0*Y**2)
     -     * ((XX+IMAG*PI)/Y+(YY+IMAG*PI)/X)
      X = T/S
      Y = U/S
      XX = DLOG(-X)
      YY = DLOG(-Y)
      IF     ( ABS(L1234(1)+L1234(2)+L1234(3)+L1234(4)).EQ.4 ) THEN
         BFDFL = HALF
      ELSEIF ( ABS(L1234(1)+L1234(2)+L1234(3)+L1234(4)).EQ.2 ) THEN
         IF ( L1234(3).EQ.L1234(4) ) THEN
            BFDFL = MPPP(X,Y,XX,YY) + MPPP(Y,X,YY,XX)
         ELSE
            BFDFL = PPMP(X,Y,XX,YY) + PPMP(Y,X,YY,XX)
         ENDIF
      ELSEIF ( ABS(L1234(1)+L1234(2)-L1234(3)-L1234(4)).EQ.4 ) THEN
         BFDFL = MMPP(X,Y,XX,YY) + MMPP(Y,X,YY,XX)
      ELSEIF ( ABS(L1234(1)-L1234(2)+L1234(3)-L1234(4)).EQ.4 ) THEN
         BFDFL = MPMP(X,Y,XX,YY)
      ELSEIF ( ABS(L1234(1)-L1234(2)-L1234(3)+L1234(4)).EQ.4 ) THEN
         BFDFL = MPMP(Y,X,YY,XX)
      ELSE
         PRINT *, "Error: BFDFL"
         STOP
      ENDIF
      RETURN
      END
C
      DOUBLE COMPLEX FUNCTION BFDFSL (S,T,U,L1234)
      IMPLICIT NONE
      DOUBLE PRECISION S,T,U
      INTEGER L1234(4)
      DOUBLE PRECISION PI,ZETA3,ZETA4,ONE,HALF,TWO,EPS
      DOUBLE COMPLEX UNIT,IMAG
      PARAMETER ( EPS = 1D-5 , PI = 3.141592654D0 )
      PARAMETER ( ZETA3 = 1.2020569D0 , ZETA4 = PI**2/90D0 )
      PARAMETER ( ONE = 1D0, HALF = 0.5D0, TWO = 2D0 )
      PARAMETER ( UNIT = DCMPLX(1D0,0D0), IMAG = DCMPLX(0D0,1D0) )
      DOUBLE PRECISION X,Y,XX,YY
      DOUBLE COMPLEX POLYLOG
      EXTERNAL POLYLOG
      DOUBLE COMPLEX MPPP,MMPP,MPMP
      MPPP(X,Y,XX,YY) = UNIT/8D0 * ( (ONE+X**2)/Y**2*((XX+IMAG*PI)**2)
     -     + HALF*(X**2+Y**2)*((XX-YY)**2+PI**2)
     -     - 4D0*(ONE/Y-X)*(XX+IMAG*PI) )
      MMPP(X,Y,XX,YY) = -TWO*X**2 * ( POLYLOG(4,-X) + POLYLOG(4,-Y)
     -     - (XX+IMAG*PI)*(POLYLOG(3,-X)+POLYLOG(3,-Y))
     -     + XX**4/12D0 - XX**3*Y/3D0 + PI**2/12D0*X*Y - 4D0/90D0*PI**4
     -     + IMAG*PI/6D0*X*(XX**2-3D0*X*Y+PI**2) )
     -     - (X-Y)*(POLYLOG(4,-X/Y)-PI**2/6D0*POLYLOG(2,-X))
     -     - X * ( 2D0*POLYLOG(3,-X) - POLYLOG(3,-X/Y) - 3D0*ZETA3
     -     - TWO*(XX+IMAG*PI)*POLYLOG(2,-X)
     -     + (XX-YY)*(POLYLOG(2,-X/Y)+XX**2)
     -     + ONE/12D0*(5D0*(XX-YY)+18D0*IMAG*PI)*((XX-YY)**2+PI**2)
     -     - TWO/3D0*XX*(XX**2+PI**2) - IMAG*PI*(YY**2+PI**2) )
     -     + (ONE-TWO*X**2)/(4D0*Y**2)*((XX+IMAG*PI)**2+PI**2)
     -     - ONE/8D0*(TWO*X*Y+3D0)*((XX-YY)**2+PI**2) + PI**2/12D0
     -     + (HALF/Y+X)*(XX+IMAG*PI) - HALF*HALF
      MPMP(X,Y,XX,YY) = -TWO*(ONE+X**2)/Y**2 * (
     -     POLYLOG(4,-X/Y) - POLYLOG(4,-Y)
     -     + HALF*(XX-TWO*YY-IMAG*PI)*(POLYLOG(3,-X)-ZETA3)
     -     + ONE/24D0*(XX**4+TWO*IMAG*PI*XX**3-4D0*XX*YY**3
     -     + YY**4 + TWO*PI**2*YY**2) + 7D0/360D0*PI**4 )
     -     - TWO*(X-ONE)/Y * ( POLYLOG(4,-X) - ZETA4
     -     - HALF*(XX+IMAG*PI)*(POLYLOG(3,-X)-ZETA3)
     -     + PI**2/6D0*(POLYLOG(2,-X)-PI**2/6D0-HALF*XX**2)
     -     - XX**4/48D0 )
     -     + (TWO*X/Y-ONE) * (POLYLOG(3,-X) - (XX+IMAG*PI)*POLYLOG(2,-X)
     -     + ZETA3 - XX**3/6D0 - PI**2/3D0*(XX+YY) )
     -     + TWO*(TWO*X/Y+ONE) * ( POLYLOG(3,-Y)
     -     + (YY+IMAG*PI)*POLYLOG(2,-X) - ZETA3+XX/4D0*(TWO*YY**2+XX**2)
     -     - XX**2/8D0*(XX+3D0*IMAG*PI) )
     -     - HALF**2*(TWO*X**2-Y**2)*((XX-YY)**2+PI**2)
     -     - HALF**2*(3D0+TWO*X/Y**2)*((XX+IMAG*PI)**2+PI**2)
     -     - (TWO-Y**2)/(4D0*X**2)*((YY+IMAG*PI)**2+PI**2) + PI**2/6D0
     -     + HALF*(TWO*X+Y**2)*((XX+IMAG*PI)/X+(YY+IMAG*PI)/X) - HALF
      X = T/S
      Y = U/S
      XX = DLOG(-X)
      YY = DLOG(-Y)
      IF     ( ABS(L1234(1)+L1234(2)+L1234(3)+L1234(4)).EQ.4 ) THEN
         BFDFSL = -3D0*HALF
      ELSEIF ( ABS(L1234(1)+L1234(2)+L1234(3)+L1234(4)).EQ.2 ) THEN
         BFDFSL = MPPP(X,Y,XX,YY) + MPPP(Y,X,YY,XX)
      ELSEIF ( ABS(L1234(1)+L1234(2)-L1234(3)-L1234(4)).EQ.4 ) THEN
         BFDFSL = MMPP(X,Y,XX,YY) + MMPP(Y,X,YY,XX)
      ELSEIF ( ABS(L1234(1)-L1234(2)+L1234(3)-L1234(4)).EQ.4 ) THEN
         BFDFSL = MPMP(X,Y,XX,YY)
      ELSEIF ( ABS(L1234(1)-L1234(2)-L1234(3)+L1234(4)).EQ.4 ) THEN
         BFDFSL = MPMP(Y,X,YY,XX)
      ELSE
         PRINT *, "Error: BFDFSL"
         STOP
      ENDIF
      RETURN
      END
C
      DOUBLE COMPLEX FUNCTION POLYLOG(N1,X)
      IMPLICIT NONE
      DOUBLE PRECISION X
      INCLUDE 'CHAPLIN.INC'
      Z = DCMPLX(X)
      IF (N1.EQ.2) THEN
         POLYLOG = HPL1(0,1,Z)
      ELSEIF (N1.EQ.2) THEN
         POLYLOG = HPL1(0,0,1,Z)
      ELSEIF (N1.EQ.2) THEN
         POLYLOG = HPL1(0,0,0,1,Z)
      ELSE
         WRITE(6,*) 'ERROR:',N1,Z
         STOP
      ENDIF
      RETURN
      END
