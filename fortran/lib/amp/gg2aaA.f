C     gg2aaA.f: the top-loop amplitude at the threshold point
C     s = 4 m_t**2 (r = 1) in closed form, as a function of cos(theta*)
C     only (independent of m_t).  Used by the 2016 partial-wave driver
C     legacy/Pwd_gg2aa.f (and declared in legacy/Argand_gg2aa.f); the
C     coefficients A0PPPP, A2PPPP, ... in gg2aaG.f are presumably its
C     partial-wave projections (inference, not recomputed here).
C     Not called by mktemplate.
C
C     GG2AAA (COS, L1234)
C       COS    cos(theta*) in the gg rest frame
C       L1234  helicities (+-1) of the four bosons, classes as in GG2AAQ
C     Returns the real part (dimensionless) of GG2AAT at r = 1, with
C       S0 = t/(4 m_t**2) = -(1-COS)/2,  T0 = u/(4 m_t**2) = -(1+COS)/2.
C     BSFN, BTFN = Re B(S0), B(T0) (bubble, BS and BT are the beta of
C     S0 and T0); TSFN, TTFN = T(S0), T(T0) (triangle); -pi**2/4 is
C     Re T(1); IRS, IRT, IST = Re I(1,S0), Re I(1,T0), Re I(S0,T0)
C     (box).  Same helicity classes and normalisation as GG2AAT.
C     DDILOG is declared but not used.
      DOUBLE PRECISION FUNCTION GG2AAA (COS,L1234)
      IMPLICIT NONE
      DOUBLE PRECISION COS
      INTEGER L1234(4)
      DOUBLE PRECISION ONE,TWO,FOUR,HALF,PI
      PARAMETER ( ONE = 1D0, TWO = 2D0, FOUR = 4D0, HALF = 0.5D0 )
      PARAMETER ( PI = 3.141592654D0 )
      DOUBLE PRECISION S0,T0,A1,A2,A3,BS,BT
      DOUBLE PRECISION BSFN,BTFN,TSFN,TTFN,FRFN,IRS,IRT,IST
      DOUBLE COMPLEX   FFN
      S0 = -(ONE-COS)*HALF
      T0 = -(ONE+COS)*HALF
      A1 = DSQRT(ONE/(-S0))
      A2 = DSQRT(ONE/(-T0))
      A3 = DSQRT((5D0-COS**2)/(ONE-COS**2))
      BS = DSQRT((3D0-COS)/(ONE-COS))
      BT = DSQRT((3D0+COS)/(ONE+COS))
      BSFN = - ONE + HALF*BS*DLOG((BS+ONE)/(BS-ONE))
      BTFN = - ONE + HALF*BT*DLOG((BT+ONE)/(BT-ONE))
      TSFN = DLOG((BS+ONE)/(BS-ONE))**2/FOUR
      TTFN = DLOG((BT+ONE)/(BT-ONE))**2/FOUR
      IRS = FRFN(A1) + DREAL(FFN(S0,BS,A1))
      IRT = FRFN(A2) + DREAL(FFN(T0,BT,A2))
      IST = DREAL(FFN(S0,BS,A3) + FFN(T0,BT,A3))
      IF     ( ABS(L1234(1)+L1234(2)+L1234(3)+L1234(4)).EQ.4 ) THEN
         GG2AAA = ONE - HALF/S0 * IRS - HALF/T0 * IRT
     -        - HALF/(S0*T0) * IST
      ELSEIF ( ABS(L1234(1)+L1234(2)+L1234(3)+L1234(4)).EQ.2 ) THEN
         GG2AAA = ONE + (ONE+ONE/S0+ONE/T0) * (-PI**2/FOUR + TSFN + TTFN)
     -        - (ONE/T0+HALF/S0)*IRS - (ONE/S0+HALF/T0)*IRT
     -        - (ONE+HALF/(S0*T0))*IST
      ELSEIF ( ABS(L1234(1)+L1234(2)-L1234(3)-L1234(4)).EQ.4 ) THEN
         GG2AAA = - ONE - TWO*(ONE+TWO*S0)*BSFN - TWO*(ONE+TWO*T0)*BTFN
     -        + FOUR*S0*T0*(TSFN+TTFN) + HALF/S0*IRS + HALF/T0*IRT
     -        + (HALF/(S0*T0)-TWO-FOUR*S0*T0)*IST
      ELSEIF ( ABS(L1234(1)-L1234(2)-L1234(3)+L1234(4)).EQ.4 ) THEN
         GG2AAA = - ONE + TWO*(ONE-S0)/T0*(ONE+BSFN)
     -        - TWO*((ONE+S0**2)/T0**2-ONE/T0)*(-PI**2/FOUR+TSFN)
     -        + (ONE-HALF/T0)*(IRT+IST/S0)
      ELSEIF ( ABS(L1234(1)-L1234(2)+L1234(3)-L1234(4)).EQ.4 ) THEN
         GG2AAA = - ONE + TWO*(ONE-T0)/S0*(ONE+BTFN)
     -        - TWO*((ONE+T0**2)/S0**2-ONE/S0)*(-PI**2/FOUR+TTFN)
     -        + (ONE-HALF/S0)*(IRS+IST/T0)
      ELSE
         STOP
      ENDIF
      RETURN
      END

C
C     FRFN (A): the R = 1 term (b = 0) of the box function I(1,S) of
C     GG2AAT/IFN in closed form, with A = sqrt(1 - (1+S)/S); real.
      DOUBLE PRECISION FUNCTION FRFN (A)
      IMPLICIT NONE
      DOUBLE PRECISION A
      DOUBLE PRECISION ONE,TWO,PI
      PARAMETER ( ONE = 1D0, TWO = 2D0 )
      PARAMETER ( PI = 3.141592654D0 )
      DOUBLE COMPLEX DL1,DL2
      INCLUDE 'CHAPLIN.INC'
      Z   = DCMPLX(A/(A+ONE),0D0)
      DL1 = HPL2(0,1,Z)
      Z   = DCMPLX((A-ONE)/A,0D0)
      DL2 = HPL2(0,1,Z)
      FRFN = ( DLOG(A)**2 + DLOG(ONE+A)**2 - TWO*DLOG(A)*DLOG(A-ONE)
     -     - TWO*PI**2/3D0 + TWO*DREAL(DL1) + TWO*DREAL(DL2) ) / (TWO*A)
      RETURN
      END
      
C     FFN (S, B, A): the S term of the box function I(1,S) (or of
C     I(S0,T0)), written with dilogarithms of (A+-1)/(A+-B) and a log
C     of (S - i EPS)(A**2 - B**2); complex, the caller keeps the real
C     part.  B is the beta of S, A as in FRFN (or A3 for I(S0,T0)).
      DOUBLE COMPLEX FUNCTION FFN (S,B,A)
      IMPLICIT NONE
      DOUBLE PRECISION S,A,B
      DOUBLE PRECISION ONE,TWO,EPS
      PARAMETER ( ONE = 1D0, TWO = 2D0, EPS = 1D-10 )
      DOUBLE COMPLEX DL1,DL2,DL3,DL4
      INCLUDE 'CHAPLIN.INC'
      Z = DCMPLX((A+ONE)/(A+B),0D0)
      DL1 = HPL2(0,1,Z)
      Z = DCMPLX((A-ONE)/(A+B),0D0)
      DL2 = HPL2(0,1,Z)
      Z = DCMPLX((A+ONE)/(A-B),0D0)
      DL3 = HPL2(0,1,Z)
      Z = DCMPLX((A-ONE)/(A-B),0D0)
      DL4 = HPL2(0,1,Z)
      FFN = ( LOG((S-DCMPLX(0D0,EPS))*(A**2-B**2)) * DLOG((A+ONE)/(A-ONE))
     -     - DL1 + DL2 - DL3 + DL4 ) / (TWO*A)
      RETURN
      END
