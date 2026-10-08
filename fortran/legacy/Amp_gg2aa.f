      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER I,J,K,L
      DOUBLE PRECISION MT,GT,MU,ASG,RS,COS
      INTEGER L1234L(4,16)
      DOUBLE COMPLEX AMPQ,AMPT,AMPG,GG2AAQ,GG2AAT,GG2AAG
      EXTERNAL                      GG2AAQ,GG2AAT,GG2AAG
      DOUBLE PRECISION Q2Q,Q2T
      DOUBLE PRECISION R,S
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
      include 'qcdparam.inc'
      include 'qcdfunc.inc'
      NQCD = 0
      CALL QCDLIBS
C.....
      Q2Q = 11D0/9D0
      Q2T =  4D0/9D0
      RS = 343D0
      MT = 173D0
      GT = 1.498D0
      MU  = 160D0
      ASG = ASQCD (MU)
      PRINT *, MU,ASG
      COS = 0.5D0
      DO I = 1,51
         COS = 0.04D0*(I-1) - 1D0
         AMPQ = GG2AAQ (   RS,COS,L1234L(1,1))
         AMPT = GG2AAT (MT,RS,COS,L1234L(1,1))
         AMPG = GG2AAG (MT,GT,MU,ASG,RS,COS,L1234L(1,1),1)
         WRITE(6,'(F10.3,4(1PE15.6))') COS, AMPT,AMPG
         WRITE(1,'(F10.3,4(1PE15.6))') COS, AMPT,AMPG
      ENDDO
      COS = 0D0
      DO I = 1,101
         RS = 1D0*(I-1) + 300D0
c     AMPQ = GG2AAQ (   RS,COS,L1234L(1,7))
         AMPT = GG2AAT (MT,RS,COS,L1234L(1,1))
         AMPG = GG2AAG (MT,GT,MU,ASG,RS,COS,L1234L(1,1),1)
         WRITE(6,'(F10.3,4(1PE15.6))') RS, AMPT,AMPG
         WRITE(2,'(F10.3,4(1PE15.6))') RS, AMPT,AMPG
         WRITE(3,'(F10.3,4(1PE15.6))') RS, AMPT+AMPG
      ENDDO
      STOP
      END
