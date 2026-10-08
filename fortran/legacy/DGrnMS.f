      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER I,J,K
      DOUBLE PRECISION MT,GT,MU,ASG
      DOUBLE PRECISION MAA,E,REG,IMG,REG0,IMG0
      DOUBLE PRECISION RR(-2:2),II(-2:2),H,DR,DI
      include 'qcdparam.inc'
      include 'qcdfunc.inc'
      NQCD = 1
      CALL QCDLIBS
C.....
      MT = 173D0
      GT = 1.5D0
C.....
      MU  = 40D0
      ASG = ASQCD (MU)
C.....
      H = 0.1D0
      DO I = -1,303
         RR(-2) = RR(-1)
         RR(-1) = RR( 0)
         RR( 0) = RR( 1)
         RR( 1) = RR( 2)
         II(-2) = II(-1)
         II(-1) = II( 0)
         II( 0) = II( 1)
         II( 1) = II( 2)
C...  
         MAA = 330D0 + H*(I-1)
         E = MAA - 2D0*MT
c     0:G0eps, 1:G1eps, 2:G0Gt, 3:G1Gt, 4:GLO
c         CALL GRNLOMSB  (E,MT,GT,4,ASG,   MU,REG,IMG)
         CALL GRNNLOMSB (E,MT,GT,ASG,MU,MU,1,REG,IMG) ! 0:LO,1:NLO
         RR( 2)  = REG
         II( 2)  = IMG
C.....
         IF ( I.GE.3 ) THEN
            DR = ( RR(-2) - 8D0*RR(-1) + 8D0*RR( 1) - RR( 2) ) / (12D0*H)
            DI = ( II(-2) - 8D0*II(-1) + 8D0*II( 1) - II( 2) ) / (12D0*H)
            PRINT *,MAA-2D0*H,REG,IMG,DR,DI
            WRITE (1,*) MAA, REG,IMG
            WRITE (2,*) MAA-2D0*H, DR,DI,(MAA-2D0*H)*DR/RR( 0),
     -           (MAA-2D0*H)*DI/II( 0)
         ELSEIF ( I.GE.1 ) THEN
            PRINT *,MAA-2D0*H,REG,IMG
            WRITE (1,*) MAA, REG,IMG
         ENDIF
      ENDDO
C.....
      STOP
      END
