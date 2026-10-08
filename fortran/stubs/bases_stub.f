C     No-op replacements for the BASES histogram/set-up routines.
C     The integrands of the gg2aa programs call XHINIT/XHFILL for
C     diagnostic histograms only; the cross sections themselves are
C     computed with VEGAS (intvegas.f), so BASES is not linked.
      SUBROUTINE XHINIT (ID,XLOW,XHIGH,NBIN,TITLE)
      IMPLICIT NONE
      INTEGER ID,NBIN
      DOUBLE PRECISION XLOW,XHIGH
      CHARACTER*(*) TITLE
      RETURN
      END
C
      SUBROUTINE XHFILL (ID,X,FX)
      IMPLICIT NONE
      INTEGER ID
      DOUBLE PRECISION X,FX
      RETURN
      END
C
C     BASES set-up calls made from the (unused) USERIN routines.
      SUBROUTINE BSSETD (NDIM,NWILD,XL,XU,IG)
      IMPLICIT NONE
      INTEGER NDIM,NWILD,IG(*)
      DOUBLE PRECISION XL(*),XU(*)
      RETURN
      END
C
      SUBROUTINE BSSETP (NCALL,ITMX1,ITMX2,ACC1,ACC2)
      IMPLICIT NONE
      INTEGER NCALL,ITMX1,ITMX2
      DOUBLE PRECISION ACC1,ACC2
      RETURN
      END
