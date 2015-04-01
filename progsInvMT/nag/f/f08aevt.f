      SUBROUTINE F08AEV(N,ALPHA,X,INCX,TAU)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DLARFG(N,ALPHA,X,INCX,TAU)
C
C  Purpose
C  =======
C
C  DLARFG generates a real elementary reflector H of order n, such
C  that
C
C        H * ( alpha ) = ( beta ),   H' * H = I.
C            (   x   )   (   0  )
C
C  where alpha and beta are scalars, and x is an (n-1)-element real
C  vector. H is represented in the form
C
C        H = I - tau * ( 1 ) * ( 1 v' ) ,
C                      ( v )
C
C  where tau is a real scalar and v is a real (n-1)-element
C  vector.
C
C  If the elements of x are all zero, then tau = 0 and H is taken to be
C  the unit matrix.
C
C  Otherwise  1 <= tau <= 2.
C
C  Arguments
C  =========
C
C  N       (input) INTEGER
C          The order of the elementary reflector.
C
C  ALPHA   (input/output) DOUBLE PRECISION
C          On entry, the value alpha.
C          On exit, it is overwritten with the value beta.
C
C  X       (input/output) DOUBLE PRECISION array, dimension
C                         (1+(N-2)*abs(INCX))
C          On entry, the vector x.
C          On exit, it is overwritten with the vector v.
C
C  INCX    (input) INTEGER
C          The increment between elements of X. INCX <> 0.
C
C  TAU     (output) DOUBLE PRECISION
C          The value tau.
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA, TAU
      INTEGER           INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  BETA, RSAFMN, SAFMIN, XNORM
      INTEGER           J, KNT
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, F06BNF, X02AMF
      EXTERNAL          DNRM2, F06BNF, X02AMF
C     .. External Subroutines ..
      EXTERNAL          DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN
C     .. Executable Statements ..
C
      IF (N.LE.1) THEN
         TAU = ZERO
         RETURN
      END IF
C
      XNORM = DNRM2(N-1,X,INCX)
C
      IF (XNORM.EQ.ZERO) THEN
C
C        H  =  I
C
         TAU = ZERO
      ELSE
C
C        general case
C
         BETA = -SIGN(F06BNF(ALPHA,XNORM),ALPHA)
         SAFMIN = X02AMF()
         IF (ABS(BETA).LT.SAFMIN) THEN
C
C           XNORM, BETA may be inaccurate; scale X and recompute them
C
            RSAFMN = ONE/SAFMIN
            KNT = 0
   20       CONTINUE
            KNT = KNT + 1
            CALL DSCAL(N-1,RSAFMN,X,INCX)
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF (ABS(BETA).LT.SAFMIN) GO TO 20
C
C           New BETA is at most 1, at least SAFMIN
C
            XNORM = DNRM2(N-1,X,INCX)
            BETA = -SIGN(F06BNF(ALPHA,XNORM),ALPHA)
            TAU = (BETA-ALPHA)/BETA
            CALL DSCAL(N-1,ONE/(ALPHA-BETA),X,INCX)
C
C           If ALPHA is subnormal, it may lose relative accuracy
C
            ALPHA = BETA
            DO 40 J = 1, KNT
               ALPHA = ALPHA*SAFMIN
   40       CONTINUE
         ELSE
            TAU = (BETA-ALPHA)/BETA
            CALL DSCAL(N-1,ONE/(ALPHA-BETA),X,INCX)
            ALPHA = BETA
         END IF
      END IF
C
      RETURN
C
C     End of F08AEV (DLARFG)
C
      END
