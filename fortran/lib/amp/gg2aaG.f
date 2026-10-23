C     gg2aaG.f: threshold part of the top-loop gg -> gamma gamma
C     amplitude, matched to the NRQCD Green function of t tbar
C     (arXiv:1607.00990):
C       M_t(match) = M_t(1-loop) + B [ G(0,0;E+i Gamma_t) - G0(E) ].
C     Near threshold the top loop is written as A + B G: A is the part
C     without the t tbar threshold singularity (a polynomial in cos via
C     partial waves), B G the t tbar Green-function part.  Only the
C     classes ++++ and --++ (helicity difference 0 in both pairs, which
C     contain J = 0, the colour-singlet 1S0 t tbar state) get a Green-
C     function term.
C
C     GG2AAG (MT, GT, MU, ASG, RS, COS, L1234, J)
C       MT, GT  top-quark mass and width [GeV]
C       MU      scale of alpha_s in the Green function [GeV]
C       ASG     alpha_s(MU)
C       RS      M_aa = sqrt(s_hat) [GeV]
C       COS     cos(theta*) in the gg rest frame
C       L1234   helicities (+-1), classes as in GG2AAQ / GG2AAT
C       J       what is returned (= MODE in mktemplate.f):
C               0: 0;  1: A + B G0;  2: A + B G;
C               3: B (G - G0), to be added to the full one-loop GG2AAT;
C               4: B G with A switched off (see the note at J.EQ.4)
C     Returns the dimensionless reduced amplitude, normalised as GG2AAT
C     (the caller multiplies by Q2T = 4/9).
C     Energy variable: E = M_aa - 2 m_t [GeV]; Gamma_t enters as
C     E + i Gamma_t inside GRNNLOMSB.
C     J=0: 0, 1: A + B*G_0, 2: A + B*G, 3: B*(G-G_0)
      DOUBLE COMPLEX FUNCTION GG2AAG (MT,GT,MU,ASG,RS,COS,L1234,J)
      IMPLICIT NONE
      DOUBLE PRECISION MT,GT,MU,ASG,RS,COS
      INTEGER L1234(4),J
      DOUBLE PRECISION ONE,TWO,FOUR,HALF
      DOUBLE PRECISION CF,PI,EPS
      DOUBLE COMPLEX IMAG
      PARAMETER ( ONE = 1D0, TWO = 2D0, FOUR = 4D0, HALF = 0.5D0 )
      PARAMETER ( CF = 1.333333333D0, PI = 3.141592654D0, EPS = 1D-12 )
      PARAMETER ( IMAG = DCMPLX(0D0,1D0) )
      INTEGER IA,IB
      DOUBLE PRECISION E,E0,R
      DOUBLE COMPLEX G0,GRN,GRN0,GRN1
C     COMMON /GRN/ caches the last Green function: E0 = its energy,
C     GRN = G, G - G0 or G0 (by J), GRN0 = G0, IA = 1 (include A) or 0
C     (drop A); IB is not used.
      COMMON /GRN/E0,GRN,GRN0,IA,IB
      DATA E0/12345.67890D0/
C     Order of the QCD potential in the Green function (0:LO, 1:NLO);
C     set through COMMON /GRNORD/, default 1 (BLOCK DATA GRNORDBD).
      INTEGER IORDG
      COMMON /GRNORD/ IORDG
      DOUBLE PRECISION ATPPPP,ATMPPP,ATMMPP,ATMPPM
      DOUBLE PRECISION A0PPPP,A2PPPP,A4PPPP,A0MMPP,A2MMPP,A4MMPP
      DOUBLE PRECISION A2MPPP,A4MPPP,A2MPPM,A3MPPM,A4MPPM
      DOUBLE PRECISION BPPPP,BMMPP
      DOUBLE PRECISION D200,D400,D220,D420,D222,D322,D422
      DOUBLE PRECISION REG,IMG
      PARAMETER ( A0PPPP=-1.0663565D0, A2PPPP=-4.97776D-3,A4PPPP=-5.3889D-5 )
      PARAMETER ( A0MMPP= 1.5738019D0, A2MMPP= 2.13711D-3,A4MMPP= 1.6898D-5 )
      PARAMETER ( A2MPPP=-2.90941D-3,A4MPPP=-1.0420D-5                   )
      PARAMETER ( A2MPPM= 1.1920027D-1,A3MPPM=-6.0737D-4,A4MPPM= 2.9467D-4 )
      DOUBLE PRECISION AT1PPPP,AT1MMPP,AT10
      INCLUDE 'qcdfunc.inc'
C     Wigner functions d^J_{lambda,lambda'}(cos) used in the partial-
C     wave expansion of A:  D200 = d^2_00, D400 = d^4_00, D220 = d^2_20,
C     D420 = d^4_20, D222 = d^2_22, D322 = d^3_22, D422 = d^4_22.
      D200 (COS) = HALF*(3D0*COS**2-ONE)
      D400 (COS) = (35D0*COS**4-30D0*COS**2+3D0)/8D0
      D220 (COS) = DSQRT(3D0/8D0)*(ONE-COS**2)
      D420 (COS) = DSQRT(5D0/32D0)*(ONE-COS**2)*(7D0*COS**2-ONE)
      D222 (COS) = (HALF*(ONE+COS))**2
      D322 (COS) = (HALF*(ONE+COS))**2*(3D0*COS-TWO)
      D422 (COS) = (HALF*(ONE+COS))**2*(7D0*COS**2-7D0*COS+ONE)
C     Free Green function without QCD and width [GeV**2],
C     R = s/(4 m_t**2):
C       G0 = -(m_t**2/(4 pi)) sqrt(1 - R - i EPS)
C     For M_aa = 2 m_t + E, m_t**2 (1 - R) = -m_t E - E**2/4, so near
C     threshold G0 = -(m_t/(4 pi)) sqrt(-m_t E).
      G0(R,EPS) = - MT**2/(FOUR*PI) * SQRT( ONE - R - IMAG*EPS )
C     A part: A(cos) = sum_J (2J+1) A_J d^J(cos), with the threshold
C     (r = 1) partial-wave coefficients A_J of the PARAMETERs above:
C     ++++ and --++ have J = 0, 2, 4 (d^J_00), -+++ has J = 2, 4
C     (d^J_20), -++- has J = 2, 3, 4 (d^J_22); -+-+ is ATMPPM(-COS).
C     A depends neither on M_aa nor on m_t.
      ATPPPP (COS) = A0PPPP + 5D0*A2PPPP*D200(COS) + 9D0*A4PPPP*D400(COS)
      ATMPPP (COS) = 5D0*A2MPPP*D220(COS) + 9D0*A4MPPP*D420(COS)
      ATMMPP (COS) = A0MMPP + 5D0*A2MMPP*D200(COS) + 9D0*A4MMPP*D400(COS)
      ATMPPM (COS) = 5D0*A2MPPM*D222(COS) + 7D0*A3MPPM*D322(COS)
     -     + 9D0*A4MPPM*D422(COS)
C     Matching coefficients B [GeV**-2] of the Green function for the
C     classes ++++ and --++; B G is dimensionless.
      BPPPP = -4D0*PI**2/MT**2
      BMMPP =  4D0*PI**2/MT**2
C     AT1PPPP, AT1MMPP: O(alpha_s) hard-coefficient terms with
C     log(MU/(2 m_t)) and the unknown constant AT10 = 0.  They are added
C     only in the commented-out lines below, so they do not enter.
C     Note (docs/REVIEW.md P17).
      AT10  = 0D0               ! not known yet. vary [-10:10].
      AT1PPPP = -MT**2/4D0*CF*BPPPP*DLOG(MU/(2D0*MT)) + AT10
      AT1MMPP = -MT**2/4D0*CF*BMMPP*DLOG(MU/(2D0*MT)) + AT10
      GG2AAG = 0D0
      IF ( J.EQ.0 ) RETURN
      R = RS**2/(4D0*MT**2)
      E = RS - 2D0*MT
C     The Green function depends on E only (for fixed MT, GT, MU, J), so
C     it is computed once per M_aa and reused for the other helicity
C     classes at the same M_aa.
C     Note (docs/REVIEW.md P8): the cache key is E alone; a change of
C     MT, GT, MU or J at the same E returns the previous GRN and IA.
      IF ( E.NE.E0 ) THEN
         GRN0 = G0(R,EPS)
         IF ( J.EQ.1 ) THEN
            IA = 1
            GRN = GRN0
         ELSE
c     0:G0eps, 1:G1eps, 2:G0Gt, 3:G1Gt, 4:GLO
c     CALL GRNLOMSB  (E,MT,GT,4,ASG,   MU,REG,IMG)
C     G(0,0;E+i Gamma_t) [GeV**2] from the Schroedinger equation with
C     the LO (IORDG = 0) or NLO (IORDG = 1) QCD potential and alpha_s =
C     ASG at the scale MU (lib/green/GrnMSBNLO.f).  J = 2 keeps A and G;
C     J = 3 drops A and keeps G - G0.
            CALL GRNNLOMSB (E,MT,GT,ASG,MU,MU,IORDG,REG,IMG) ! 0:LO,1:NLO
            IF ( J.EQ.2 ) THEN
               IA = 1
               GRN = DCMPLX(REG,IMG)
            ELSEIF ( J.EQ.3 ) THEN
               IA = 0
               GRN = DCMPLX(REG,IMG) - GRN0
            ENDIF
         ENDIF
         E0 = E
C     J = 4 is meant to follow a call with J = 2 or 3 at the same M_aa
C     (see the commented-out calls in legacy/DSig_gg2aa.f): on a cache
C     hit it returns B GRN with A switched off and then restores IA = 1.
C     On a cache miss J = 4 computes G but stores neither G nor IA (only
C     J = 2, 3 do), and GRN1 is not set.  Not used by mktemplate.
      ELSEIF ( J.EQ.4 ) THEN
         IA = 0
         GRN1 = GRN
c         GRN = DCMPLX(0.,DIMAG(GRN))
      ENDIF
C     Assemble the class: IA * A(COS) plus B GRN for ++++ and --++.
C     The ASG/PI*AT1 terms are commented out.
      IF     ( ABS(L1234(1)+L1234(2)+L1234(3)+L1234(4)).EQ.4 ) THEN
         GG2AAG = IA*ATPPPP ( COS) + BPPPP * GRN
c         GG2AAG = GG2AAG + ASG/PI*AT1PPPP
      ELSEIF ( ABS(L1234(1)+L1234(2)+L1234(3)+L1234(4)).EQ.2 ) THEN
         GG2AAG = IA*ATMPPP ( COS)
      ELSEIF ( ABS(L1234(1)+L1234(2)-L1234(3)-L1234(4)).EQ.4 ) THEN
         GG2AAG = IA*ATMMPP ( COS) + BMMPP * GRN
c         GG2AAG = GG2AAG + ASG/PI*AT1MMPP
      ELSEIF ( ABS(L1234(1)-L1234(2)-L1234(3)+L1234(4)).EQ.4 ) THEN
         GG2AAG = IA*ATMPPM ( COS)
      ELSEIF ( ABS(L1234(1)-L1234(2)+L1234(3)-L1234(4)).EQ.4 ) THEN
         GG2AAG = IA*ATMPPM (-COS)
      ELSE
         STOP
      ENDIF
      
      IF ( J.EQ.4 ) THEN
         GRN = GRN1
         IA  = 1
      ENDIF
      
      RETURN
      END
C
C     Default of IORDG (NLO potential); mktemplate.f sets it from the
C     namelist.
      BLOCK DATA GRNORDBD
      INTEGER IORDG
      COMMON /GRNORD/ IORDG
      DATA IORDG /1/
      END
