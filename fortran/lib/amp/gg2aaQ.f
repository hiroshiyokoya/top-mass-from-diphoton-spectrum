C     gg2aaQ.f: one-loop gg -> gamma gamma amplitude of a massless
C     quark loop (the light quarks u, d, s, c, b).
C
C     GG2AAQ (RS, COS, L1234)
C       RS     invariant mass M_aa = sqrt(s_hat) [GeV]
C       COS    cos(theta*) of the photons in the gg rest frame
C       L1234  helicities (+1/-1) of the four external bosons; entries
C              (1,2) and (3,4) form the two pairs (gluons, photons)
C     Returns the dimensionless reduced amplitude for one quark of unit
C     charge, normalised so that the ++++ and -+++ amplitudes are 1.
C     The caller multiplies by sum e_q**2 (Q2Q = 11/9 in mktemplate.f);
C     couplings, colour and the overall normalisation are in the
C     cross-section factor there.  Mandelstam variables in the code:
C       s = RS**2,  t = -s(1-COS)/2,  u = -s(1+COS)/2.
C
C     Helicity classes (the test on L1234 picks the class; a set and
C     the set with all signs flipped belong to the same class):
C       |l1+l2+l3+l4| = 4   ++++          1
C       |l1+l2+l3+l4| = 2   -+++ (8 sets) 1
C       |l1+l2-l3-l4| = 4   --++          MQMMPP(s,t,u)
C       |l1-l2-l3+l4| = 4   -++-          MQMPPM(s,t,u)
C       |l1-l2+l3-l4| = 4   -+-+          MQMPPM(s,u,t)
C       anything else       STOP
      DOUBLE COMPLEX FUNCTION GG2AAQ (RS,COS,L1234)
      IMPLICIT NONE
      DOUBLE PRECISION RS,COS
      INTEGER L1234(4)
      DOUBLE PRECISION PI,ONE,HALF,EPS
      DOUBLE COMPLEX IMAG
      PARAMETER ( EPS = 1D-5 , PI = 3.141592654D0 )
      PARAMETER ( ONE = 1D0, HALF = 0.5D0, IMAG = DCMPLX(0D0,1D0) )
      DOUBLE PRECISION S,T,U
      DOUBLE COMPLEX MQPPPP,MQMPPP,MQMMPP,MQMPPM
      DATA MQPPPP /1D0/, MQMPPP /1D0/
C     Statement functions for the --++ and -++- classes.  MQMMPP is
C     real for t, u < 0; MQMPPM carries an explicit imaginary part
C     -i pi (...).  EPS is not used.
      MQMMPP (S,T,U) = - HALF*(T**2+U**2)/S**2 * ( DLOG(T/U)**2
     -     + PI**2 ) - (T-U)/S * DLOG(T/U) - ONE
      MQMPPM (S,T,U) = - HALF*(T**2+S**2)/U**2 * LOG(-T/S)**2
     -     - (T-S)/U * LOG(-T/S) - ONE - IMAG*PI * (
     -     (T**2+S**2)/U**2 * LOG(-T/S) + (T-S)/U )
c     IF ( COS.EQ. 1D0 ) COS =  1D0-EPS
c     IF ( COS.EQ.-1D0 ) COS = -1D0+EPS
      S = RS**2
      T = - HALF*S*(ONE-COS)
      U = - HALF*S*(ONE+COS)
      IF     ( ABS(L1234(1)+L1234(2)+L1234(3)+L1234(4)).EQ.4 ) THEN
         GG2AAQ = MQPPPP
      ELSEIF ( ABS(L1234(1)+L1234(2)+L1234(3)+L1234(4)).EQ.2 ) THEN
         GG2AAQ = MQMPPP
      ELSEIF ( ABS(L1234(1)+L1234(2)-L1234(3)-L1234(4)).EQ.4 ) THEN
         GG2AAQ = MQMMPP (S,T,U)
      ELSEIF ( ABS(L1234(1)-L1234(2)-L1234(3)+L1234(4)).EQ.4 ) THEN
         GG2AAQ = MQMPPM (S,T,U)
c     GG2AAQ = MQMMPP (U,T,S)
      ELSEIF ( ABS(L1234(1)-L1234(2)+L1234(3)-L1234(4)).EQ.4 ) THEN
         GG2AAQ = MQMPPM (S,U,T)
c     GG2AAQ = MQMMPP (U,S,T)
      ELSE
         STOP
      ENDIF
      RETURN
      END
