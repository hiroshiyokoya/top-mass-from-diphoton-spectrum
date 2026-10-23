C     GrnMSBLO.f: analytic forms of the t tbar Green function at the
C     origin, G(0,0;E), for the simple cases selected by I.  No
C     Schroedinger equation is solved here (for the LO/NLO potential see
C     GRNNLOMSB in GrnMSBNLO.f).  Not called by mktemplate: the call in
C     gg2aaG.f is commented out.
C
C     GRNLOMSB (E, M, GAM, I, AS, MU, REG, IMG)
C       E         energy from threshold, M_aa - 2 m_t [GeV]
C       M         top-quark mass m_t [GeV]
C       GAM       top width [GeV]; replaced by 1D-12*M for I = 0, 1
C       I         0: free G0, no width      1: G0 + O(alpha_s), no width
C                 2: free G0 with width     3: G0 + O(alpha_s), width
C                 4: 1S pole only (Breit-Wigner, see below)
C       AS        alpha_s;  MU  scale [GeV] in the O(alpha_s) log
C       REG, IMG  real and imaginary part of G [GeV**2] (output)
C
C     With p1 + i p2 = sqrt(m_t (E + i Gamma)):
C       G0 = (m_t/(4 pi)) (-p2 + i p1) = -(m_t/(4 pi)) sqrt(-m_t(E+iG))
C     O(alpha_s) term (one Coulomb exchange), P0 = C_F alpha_s m_t/2:
C       + (m_t P0/(4 pi)) log(m_t**2 D**2/(p1**2 + p2**2))
C       + i (m_t P0/(2 pi)) atan2(p1, p2),   D = MU/(2 m_t) exp(1/2)
C     (D is marked MSB in the code).
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
C     I = 4: only the 1S Coulomb pole (the loop runs to N = 1; the sum
C     to 300 is cut off),
C       G = -|psi_1(0)|**2/(E - E_1 + i GAM),
C       |psi_n(0)|**2 = P0**3/(pi n**3),  E_n = -P0**2/(m_t n**2).
C     REG and IMG are overwritten, so the G0 and O(alpha_s) parts above
C     are not included (docs/REVIEW.md P22).
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
