      SUBROUTINE F08MEX(F,G,H,SSMIN,SSMAX)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DLAS2(F,G,H,SSMIN,SSMAX)
C
C  Purpose
C  =======
C
C  DLAS2  computes the singular values of the 2-by-2 matrix
C     [  F   G  ]
C     [  0   H  ].
C  On return, SSMIN is the smaller singular value and SSMAX is the
C  larger singular value.
C
C  Arguments
C  =========
C
C  F       (input) DOUBLE PRECISION
C          The (1,1) entry of the 2-by-2 matrix.
C
C  G       (input) DOUBLE PRECISION
C          The (1,2) entry of the 2-by-2 matrix.
C
C  H       (input) DOUBLE PRECISION
C          The (2,2) entry of the 2-by-2 matrix.
C
C  SSMIN   (output) DOUBLE PRECISION
C          The smaller singular value.
C
C  SSMAX   (output) DOUBLE PRECISION
C          The larger singular value.
C
C  Further Details
C  ===============
C
C  Barring over/underflow, all output quantities are correct to within
C  a few units in the last place (ulps), even in the absence of a guard
C  digit in addition/subtraction.
C
C  In IEEE arithmetic, the code works correctly if one matrix entry is
C  infinite.
C
C  Overflow will not occur unless the largest singular value itself
C  overflows, or is within a few ulps of overflow. (On machines with
C  partial overflow, like the Cray, overflow may occur if the largest
C  singular value is within a factor of 2 of overflow.)
C
C  Underflow is harmless if underflow is gradual. Otherwise, results
C  may correspond to a matrix modified by perturbations of size near
C  the underflow threshold.
C
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D0)
      DOUBLE PRECISION  TWO
      PARAMETER         (TWO=2.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  F, G, H, SSMAX, SSMIN
C     .. Local Scalars ..
      DOUBLE PRECISION  AS, AT, AU, C, FA, FHMN, FHMX, GA, HA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      FA = ABS(F)
      GA = ABS(G)
      HA = ABS(H)
      FHMN = MIN(FA,HA)
      FHMX = MAX(FA,HA)
      IF (FHMN.EQ.ZERO) THEN
         SSMIN = ZERO
         IF (FHMX.EQ.ZERO) THEN
            SSMAX = ZERO
         ELSE
            SSMAX = MAX(FHMX,GA)*SQRT(ONE+(MIN(FHMX,GA)/MAX(FHMX,GA))
     *              **2)
         END IF
      ELSE
         IF (GA.LT.FHMX) THEN
            AS = ONE + FHMN/FHMX
            AT = (FHMX-FHMN)/FHMX
            AU = (GA/FHMX)**2
            C = TWO/(SQRT(AS*AS+AU)+SQRT(AT*AT+AU))
            SSMIN = FHMN*C
            SSMAX = FHMX/C
         ELSE
            AU = FHMX/GA
            IF (AU.EQ.ZERO) THEN
C
C              Avoid possible harmful underflow if exponent range
C              asymmetric (true SSMIN may not underflow even if
C              AU underflows)
C
               SSMIN = (FHMN*FHMX)/GA
               SSMAX = GA
            ELSE
               AS = ONE + FHMN/FHMX
               AT = (FHMX-FHMN)/FHMX
               C = ONE/(SQRT(ONE+(AS*AU)**2)+SQRT(ONE+(AT*AU)**2))
               SSMIN = (FHMN*C)*AU
               SSMIN = SSMIN + SSMIN
               SSMAX = GA/(C+C)
            END IF
         END IF
      END IF
      RETURN
C
C     End of F08MEX (DLAS2)
C
      END
