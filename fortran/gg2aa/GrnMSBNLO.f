      SUBROUTINE GRNNLOMSB (E,MT,GAMT,ASB,MUB,MU,IORD,REG,IMG)
C     Green Fnc. in the MSbar scheme
C     IORD:  0:LO, 1:NLO
      IMPLICIT NONE
      DOUBLE PRECISION E,MT,GAMT,ASB,MUB,MU,REG,IMG
      INTEGER IORD,ID
      INTEGER N
      PARAMETER (N=2)
      INTEGER NLOOP
      DATA NLOOP/10000/
      DOUBLE PRECISION L10,PI,GAME
      PARAMETER (L10=2.3025851D0,PI=3.14159254D0,GAME=0.57721 5665D0)
      DOUBLE COMPLEX CZERO,CONE,CIMG
      PARAMETER (CZERO=(0D0,0D0),CONE=(1D0,0D0), CIMG=(0D0,1D0))
      INTEGER I,J,K
      DOUBLE PRECISION MM,AS,M,G,ET
      COMMON/MSB1/     MM,AS,M,G,ET,ID
      INTEGER ICLR,IPOT
      DOUBLE PRECISION CF,CA,TR
      INTEGER NF
      PARAMETER (CF=1.333333D0,CA=3D0,TR=0.5D0,NF=5)
      DOUBLE PRECISION B0,A1
      DOUBLE PRECISION H,X,R,EPS,ERR
      DOUBLE PRECISION X0,X1,R0
      COMMON/X01/X0,X1
      DATA X0,X1/-8D0,0D0/
      DOUBLE PRECISION L1
      DOUBLE COMPLEX U10(N),U1(N),U1N(N),DU1
      DOUBLE COMPLEX U00(N),U0(N),U0N(N),DU0
      DOUBLE COMPLEX B,BP,CMSB
      DOUBLE PRECISION DEL,D0
      INTEGER IB
      MM = MUB
      AS = ASB
      M  = MT
      G  = GAMT
      ET = E
      ID = IORD
C.....
      EPS = 1D-3
      ERR = EPS/DBLE(NLOOP)
C.....
      B0 = 11D0/3D0*CA - 4D0/3D0*TR*NF
      A1 = 31D0/9D0*CA - 20D0/9D0*TR*NF
C.....
c      XIN = -8D0
c      XIR =  0D0
C...  Initial Conditions
      R0 = 10D0**X0
      L1 = DLOG(MU*R0)+GAME
C...  Initial condition for u0
      U00(1) =  CONE*( 1D0 - MT*CF*ASB*R0 *
     -     ( L1 + IORD*ASB/(4D0*PI)*( (A1-2D0*B0)*L1 + B0*L1**2 ) ) )
      U00(2) = -CONE*MT*CF*ASB * ( 1D0 + L1 + IORD*ASB/(4D0*PI)*(
     -     A1-2D0*B0 + 2D0*B0*L1 + (A1-2D0*B0)*L1 + B0*L1**2 ) ) * L10 * R0
C...  Initial-condition for u1
      U10(1) = DCMPLX( R0,0D0)
      U10(2) = DCMPLX(1D0,0D0) * L10 * R0
C...  
      CMSB = MT*CF*ASB* (-0.5D0) * CONE
C...  Set Final Point of Evaluation
c     X0 = XIN
c     X1 = XIR
 101  CONTINUE
C...  Set-Initialize
      H = (X1-X0)/NLOOP         ! Step of Evaluations
      IB = 0
      D0 = ERR
      X = X0
      DO I=1,N
         U0(I) = U00(I)
         U1(I) = U10(I)
      ENDDO
C...  Take R->Infinity Limit
      DO I = 1,NLOOP,1
         CALL RUNGE_KUTTA (X,H,U0,U0N)
         CALL RUNGE_KUTTA (X,H,U1,U1N)
         X = X+H
         R = 10D0**X
         B = -U0N(1)/U1N(1)
         DEL = ABS((B-BP)/BP)
         IF ( DEL.LE.D0 ) THEN
            IB = IB*2
            D0 = DEL
         ELSE
            IB = 1
            D0 = ERR
         ENDIF
         IF ( IB.GT.1000 ) THEN
c     WRITE (1,*) E, X, "; R->oo converge"
            GOTO 110
         ENDIF
         DO K=1,N
            U0(K) = U0N(K)      ! Set GG for the next step
            U1(K) = U1N(K)      ! Set GG for the next step
         ENDDO
         BP = B
      ENDDO
c     WRITE (1,*) E, X1, "; R->oo, NOT CONVERGE"
      NLOOP = NLOOP + 5000
      X1 = X1 + .2
      GOTO 101
 110  CONTINUE
      B = (B + CMSB) * MT/(4D0*PI)
      IMG = DIMAG(B)
      REG = DREAL(B)
C.....
      RETURN
      END
C
      SUBROUTINE RUNGE_KUTTA (X,H,GG,GGN)
      IMPLICIT NONE
      INTEGER N
      PARAMETER (N=2)
      INTEGER J
      DOUBLE PRECISION X,H
      DOUBLE COMPLEX GG(N),GGN(N)
      DOUBLE COMPLEX KK1(N),KK2(N),KK3(N),KK4(N)
      CALL K1FUNC (X,H,GG,    KK1)
      CALL K2FUNC (X,H,GG,KK1,KK2)
      CALL K3FUNC (X,H,GG,KK2,KK3)
      CALL K4FUNC (X,H,GG,KK3,KK4)
      DO J=1,N
         GGN(J) = GG(J) + 1D0/6D0 *
     -        (KK1(J) + 2D0*KK2(J) + 2D0*KK3(J) + KK4(J))
      ENDDO
      RETURN
      END
C
      SUBROUTINE RKFUNC (X,H,GG,FF)
      IMPLICIT NONE
      INTEGER N
      PARAMETER (N=2)
      DOUBLE PRECISION L10
      PARAMETER (L10=2.3025851D0)
      DOUBLE PRECISION X,H,R
      DOUBLE COMPLEX GG(N), FF(N)
      DOUBLE COMPLEX DFUNC
      EXTERNAL DFUNC
      R = 10D0**X
      FF(1) = H * GG(2)
      FF(2) = H * (L10*GG(2) - L10**2*R**2 * DFUNC (X,GG(1)))
      RETURN
      END
C
      DOUBLE COMPLEX FUNCTION DFUNC (X,G1)
      IMPLICIT NONE
      DOUBLE COMPLEX CONE,CIMG
      PARAMETER (CONE=(1D0,0D0), CIMG=(0D0,1D0))
      DOUBLE PRECISION X
      DOUBLE COMPLEX G1
      DOUBLE PRECISION MU,AS,MT,GAMT,E
      INTEGER ID
      COMMON/MSB1/     MU,AS,MT,GAMT,E,ID
      DOUBLE PRECISION VQCD
      EXTERNAL VQCD
      DFUNC = MT * (E + CIMG*GAMT - VQCD(X)) * G1
      RETURN
      END
C
      SUBROUTINE K1FUNC (X,H,GG,KK1)
      IMPLICIT NONE
      INTEGER N
      PARAMETER (N=2)
      INTEGER I
      DOUBLE PRECISION X,H,X1
      DOUBLE COMPLEX GG(N),GG1(N),KK1(N)
      X1 = X
      DO I=1,N
         GG1(I) = GG(I)
      ENDDO
      CALL RKFUNC (X1,H,GG1,KK1)
      RETURN
      END
C
      SUBROUTINE K2FUNC (X,H,GG,KK1,KK2)
      IMPLICIT NONE
      INTEGER N
      PARAMETER (N=2)
      INTEGER I
      DOUBLE PRECISION X,H,X2
      DOUBLE COMPLEX GG(N),GG2(N),KK1(N),KK2(N)
      X2 = X + H/2D0
      DO I=1,N
         GG2(I) = GG(I) + KK1(I)/2D0
      ENDDO
      CALL RKFUNC (X2,H,GG2,KK2)
      RETURN
      END
C
      SUBROUTINE K3FUNC (X,H,GG,KK2,KK3)
      IMPLICIT NONE
      INTEGER N
      PARAMETER (N=2)
      INTEGER I
      DOUBLE PRECISION X,H,X3
      DOUBLE COMPLEX GG(N),GG3(N),KK2(N),KK3(N)
      X3 = X + H/2D0
      DO I=1,N
         GG3(I) = GG(I) + KK2(I)/2D0
      ENDDO
      CALL RKFUNC (X3,H,GG3,KK3)
      RETURN
      END
C
      SUBROUTINE K4FUNC (X,H,GG,KK3,KK4)
      IMPLICIT NONE
      INTEGER N
      PARAMETER (N=2)
      INTEGER I
      DOUBLE PRECISION X,H,X4
      DOUBLE COMPLEX GG(N),GG4(N),KK3(N),KK4(N)
      X4 = X + H
      DO I=1,N
         GG4(I) = GG(I) + KK3(I)
      ENDDO
      CALL RKFUNC (X4,H,GG4,KK4)
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION VQCD (X)
      IMPLICIT NONE
      DOUBLE PRECISION X,R
      DOUBLE PRECISION MUB,ASB,MT,GAMT,E
      INTEGER IORD
      COMMON/MSB1/     MUB,ASB,MT,GAMT,E,IORD
      INTEGER NF
      PARAMETER (NF=5)
      DOUBLE PRECISION PI,SQ2,EGAM, C,CF,CA
      PARAMETER (PI=3.141593D0,SQ2=1.41421356D0,EGAM=0.577216D0,
     -     CF=1.33333333D0,CA=3D0)
      DOUBLE PRECISION A1,B0
      DOUBLE PRECISION GF,MH
      PARAMETER (GF=1.16637D-5,MH=120D0)
      DOUBLE PRECISION VH
      R = 10D0**X
      B0 = 11D0/3D0*CA - 2D0/3D0*NF
      A1 = 31D0/9D0*CA - 10D0/9D0*NF
      VQCD = -CF * ASB/R
     -     * (1D0 + IORD*ASB/(4D0*PI)*(2D0*B0*(DLOG(MUB*R) + EGAM) + A1) )
c     VH = -GF*MT**2/2D0/SQ2/PI * DEXP(-MH*R) / R
c     VQCD = VQCD + VH      RETURN
      RETURN
      END
