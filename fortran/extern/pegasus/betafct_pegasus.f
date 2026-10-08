*
* ..File: betafct.f
*
*
* ..The subroutine BETAFCT for the coefficients  BETA0...BETA3  of the 
*    beta function of QCD up to order alpha_s^5 (N^3LO), normalized by 
*
*        d a_s / d ln mu_r^2  =  - BETA0 a_s^2 - BETA1 a_s^3 - ... 
*
*    with  a_s = alpha_s/(4*pi). 
*
* ..The MSbar coefficients are written to the common-block  BETA  for 
*   NF = 3...6  (parameters NFMIN, NFMAX) quark flavours.
*
* ..The factors CF, CA and TF  are taken from the common-block  COLOUR. 
*    Beyond NLO the QCD colour factors are hard-wired in this routine,
*    and the numerical coefficients are truncated to six digits.
*
* =====================================================================
*
*
       SUBROUTINE BETAFCT_PEGASUS
*
       IMPLICIT DOUBLE PRECISION (A - Z)
       INTEGER NFMIN, NFMAX, NF
       PARAMETER (NFMIN = 3, NFMAX = 6)
*
* ---------------------------------------------------------------------
*
* ..Input common-block
*
       COMMON / COLOUR / CF, CA, TR
*
* ..Output common-block
*
       COMMON / BETA   / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     1                   BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)
*
* ---------------------------------------------------------------------
*
* ..The full LO and NLO coefficients 
*
       B00 =  11./3.D0 * CA
       B01 =  -4./3.D0 * TR
       B10 =  34./3.D0 * CA**2
       B11 = -20./3.D0 * CA*TR - 4.* CF*TR
*
* ..Flavour-number loop and output to the array
*
       DO 1 NF = NFMIN, NFMAX
*
       BETA0(NF) = B00 + B01 * NF
       BETA1(NF) = B10 + B11 * NF
*
       BETA2(NF) = 1428.50 - 279.611 * NF + 6.01852 * NF**2
       BETA3(NF) = 29243.0 - 6946.30 * NF + 405.089 * NF**2 
     1             + 1.49931 * NF**3
*
* ---------------------------------------------------------------------
*
  1    CONTINUE
*
       RETURN
       END
*
* =================================================================av==
