C     gg2aaT.f: one-loop gg -> gamma gamma amplitude of the top-quark
C     loop with the full m_t dependence.
C
C     GG2AAT (MT, RS, COS, L1234)
C       MT     top-quark mass [GeV]
C       RS     M_aa = sqrt(s_hat) [GeV]
C       COS    cos(theta*) in the gg rest frame
C       L1234  helicities (+-1) of the four bosons, classes as in GG2AAQ
C     Returns the dimensionless reduced amplitude for unit charge (the
C     caller multiplies by e_t**2 = Q2T = 4/9).  Normalised like GG2AAQ:
C     MTPPPP and MTMPPP tend to 1 for m_t -> 0 (inference).  It is
C     complex above r = s/(4 m_t**2) = 1, i.e. M_aa > 2 m_t, where the
C     t tbar cut opens.
C
C     The amplitudes are written in the scaled invariants
C       r = s/(4 m_t**2),  s' = t/(4 m_t**2),  t' = u/(4 m_t**2)
C     (in the code S and T are overwritten with s' and t') and in terms
C     of the one-loop functions below:
C       BFN(r)    bubble,   e.g. beta acosh(sqrt r) - 1 - i pi beta/2
C       TFN(r)    triangle, e.g. acosh(sqrt r)**2 - pi**2/4
C                                - i pi acosh(sqrt r)
C                 (both for r >= 1, beta = sqrt(1 - 1/r))
C       IFN(r,s)  box, built from dilogarithms (CHAPLIN HPL2)
C     Helicity classes: ++++ MTPPPP(r,s',t'), -+++ MTMPPP(r,s',t'),
C     --++ MTMMPP(r,s',t'), and by crossing -++- MTMMPP(t',s',r),
C     -+-+ MTMMPP(s',t',r).
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
      
C
C     BFN (R): bubble function of the scaled invariant R (s/(4 m**2)
C     or t/(4 m**2) etc.), dimensionless:
C       R >= 1     : beta acosh(sqrt R) - 1 - i pi beta/2,
C                    beta = sqrt(1 - 1/R)
C       0 < R < 1  : beta' asin(sqrt R) - 1, beta' = sqrt(1/R - 1)
C       R = 0      : 0
C       R < 0      : beta asinh(sqrt(-R)) - 1
      DOUBLE COMPLEX FUNCTION BFN (R)
      IMPLICIT NONE
      DOUBLE PRECISION R,B
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
      
C
C     TFN (R): triangle function, dimensionless:
C       R >= 1     : acosh(sqrt R)**2 - pi**2/4 - i pi acosh(sqrt R)
C       0 <= R < 1 : -asin(sqrt R)**2
C       R < 0      : asinh(sqrt(-R))**2
      DOUBLE COMPLEX FUNCTION TFN (R)
      IMPLICIT NONE
      DOUBLE PRECISION R
      DOUBLE PRECISION ZERO,ONE, PI
      PARAMETER ( ZERO = 0D0, ONE = 1D0)
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
      
C
C     IFN (R, S): box function, dimensionless.  With
C     A = sqrt(1 - (R+S)/(R S)) it is the sum of one term for R and one
C     for S, each
C       [ -Li2((A+1)/(A+b)) + Li2((A-1)/(A+b))
C         -Li2((A+1)/(A-b)) + Li2((A-1)/(A-b)) ] / (2A),
C     with b = sqrt(1 - 1/R) (R < 0 or R > 1) or i sqrt(1/R - 1)
C     (0 <= R <= 1), and likewise for S.  Li2 = CHAPLIN HPL2(0,1,z).
C     Only the real part of the dilogarithm sum is kept; for R > 1
C     (S > 1) the imaginary part i pi/(2A) log((A-b)/(A+b)) is added.
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
