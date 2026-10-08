      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER I,J,K
      DOUBLE PRECISION MT,GT,MU,ASG
      DOUBLE PRECISION E,MAA,REG,IMG,REG0,IMG0
      include 'qcdparam.inc'
      include 'qcdfunc.inc'
      NQCD = 0
      CALL QCDLIBS
C.....
      MT = 173.0D0
      GT = 1.498D0
C.....
      MU  = 40D0
      ASG = ASQCD (MU)
C.....
      DO I = 1,1001
         E = -8D0 + 0.01D0*(I-1)
         MAA = 2D0*MT + E
c     0:G0eps, 1:G1eps, 2:G0Gt, 3:G1Gt, 4:GLO
         CALL GRNLOMSB  (E,MT,GT,4,ASG,   MU,REG0,IMG0)
c         CALL GRNNLOMSB (E,MT,GT,ASG,MU,MU,1,REG,IMG) ! 0:LO,1:NLO
         PRINT *,E,IMG0,IMG
         WRITE (11,*) E, REG0,IMG0,REG,IMG
         WRITE (12,*) MAA, REG0,IMG0,REG,IMG
      ENDDO
C.....
      STOP
      END
