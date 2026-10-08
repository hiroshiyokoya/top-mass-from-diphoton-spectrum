*
* ..File: asrgkt.f      (requires a previous call of BETAFCT)
*
*
* ..The running coupling of QCD,  
*
*         AS  =  a_s  =  alpha_s(mu_r^2)/(4 pi),
*
*    obtained by integrating the evolution equation for a fixed number
*    of massless flavours  NF.  Except at leading order (LO),  AS  is 
*    obtained using a fourth-order Runge-Kutta integration. 
*
* ..The initial and final scales  R20  and  R2,  the value  AS0  at
*    R20, and  NF  are passed as function arguments.  The coefficients 
*    of the beta function up to  a_s^5 (N^3LO)  are provided by the 
*    common-block  BETA.  The order of the expansion  NAORD  (defined 
*    as the 'n' in N^nLO) and the number of steps  NASTPS  for the 
*    integration beyond LO are given by the common-block  ASPAR.
*
* =====================================================================
*
*
      DOUBLE PRECISION FUNCTION AS_PEGASUS (R2, R20, AS0, NF)
*
       IMPLICIT DOUBLE PRECISION (A - Z)
       INTEGER NFMIN, NFMAX, NF, NAORD, NASTPS, K1
       PARAMETER (NFMIN = 3, NFMAX = 6)
       PARAMETER ( SXTH = 0.16666 66666 66666 D0 )
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks 
*
       COMMON / ASPAR  / NAORD, NASTPS
       COMMON / BETA   / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     ,                   BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)
*
* ..The beta functions FBETAn at N^nLO for n = 1, 2, and 3
*
       FBETA1(A) = - A**2 * ( BETA0(NF) + A *   BETA1(NF) )
       FBETA2(A) = - A**2 * ( BETA0(NF) + A * ( BETA1(NF)
     ,                        + A * BETA2(NF) ) )
       FBETA3(A) = - A**2 * ( BETA0(NF) + A * ( BETA1(NF)
     ,                        + A * (BETA2(NF) + A * BETA3(NF)) ) )
*
* ---------------------------------------------------------------------
*
* ..Initial value, evolution distance and step size
*
       AS = AS0
       LRRAT = LOG (R2/R20)
       DLR = LRRAT / NASTPS
*
* ..Solution of the evolution equation depending on  NAORD
*   (fourth-order Runge-Kutta beyond the leading order)
*
       IF (NAORD .EQ. 0) THEN
*
         AS = AS0 / (1.+ BETA0(NF) * AS0 * LRRAT)
*
       ELSE IF (NAORD .EQ. 1) THEN
*
       DO 2 K1 = 1, NASTPS
         XK0 = DLR * FBETA1 (AS)
         XK1 = DLR * FBETA1 (AS + 0.5 * XK0)
         XK2 = DLR * FBETA1 (AS + 0.5 * XK1)
         XK3 = DLR * FBETA1 (AS + XK2)
         AS = AS + SXTH * (XK0 + 2.* XK1 + 2.* XK2 + XK3)
  2    CONTINUE
*
       ELSE IF (NAORD .EQ. 2) THEN
*
       DO 3 K1 = 1, NASTPS
         XK0 = DLR * FBETA2 (AS)
         XK1 = DLR * FBETA2 (AS + 0.5 * XK0)
         XK2 = DLR * FBETA2 (AS + 0.5 * XK1)
         XK3 = DLR * FBETA2 (AS + XK2)
         AS = AS + SXTH * (XK0 + 2.* XK1 + 2.* XK2 + XK3)
  3    CONTINUE
*  
       ELSE IF (NAORD .EQ. 3) THEN
*
       DO 4 K1 = 1, NASTPS
         XK0 = DLR * FBETA3 (AS)
         XK1 = DLR * FBETA3 (AS + 0.5 * XK0)
         XK2 = DLR * FBETA3 (AS + 0.5 * XK1)
         XK3 = DLR * FBETA3 (AS + XK2)
         AS = AS + SXTH * (XK0 + 2.* XK1 + 2.* XK2 + XK3)
  4    CONTINUE
       END IF
*
* ---------------------------------------------------------------------
       AS_PEGASUS = AS
* ---------------------------------------------------------------------
*
       RETURN
       END
*
* =================================================================av==
