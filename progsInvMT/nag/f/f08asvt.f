      SUBROUTINE F08ASV(N,ALPHA,X,INCX,TAU)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZLARFG(N,ALPHA,X,INCX,TAU)
C
C  Purpose
C  =======
C
C  ZLARFG generates a complex elementary reflector H of order n, such
C  that
C
C        H' * ( alpha ) = ( beta ),   H' * H = I.
C             (   x   )   (   0  )
C
C  where alpha and beta are scalars, with beta real, and x is an
C  (n-1)-element complex vector. H is represented in the form
C
C        H = I - tau * ( 1 ) * ( 1 v' ) ,
C                      ( v )
C
C  where tau is a complex scalar and v is a complex (n-1)-element
C  vector. Note that H is not hermitian.
C
C  If the elements of x are all zero and alpha is real, then tau = 0
C  and H is taken to be the unit matrix.
C
C  Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .
C
C  Arguments
C  =========
C
C  N       (input) INTEGER
C          The order of the elementary reflector.
C
C  ALPHA   (input/output) COMPLEX*16
C          On entry, the value alpha.
C          On exit, it is overwritten with the value beta.
C
C  X       (input/output) COMPLEX*16 array, dimension
C                         (1+(N-2)*abs(INCX))
C          On entry, the vector x.
C          On exit, it is overwritten with the vector v.
C
C  INCX    (input) INTEGER
C          The increment between elements of X. INCX <> 0.
C
C  TAU     (output) COMPLEX*16
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
      COMPLEX*16        ALPHA, TAU
      INTEGER           INCX, N
C     .. Array Arguments ..
      COMPLEX*16        X(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
      INTEGER           J, KNT
      LOGICAL           DIVFLG
C     .. External Functions ..
      COMPLEX*16        F06CLF
      DOUBLE PRECISION  DZNRM2, F08ASU, X02AMF
      EXTERNAL          F06CLF, DZNRM2, F08ASU, X02AMF
C     .. External Subroutines ..
      EXTERNAL          ZDSCAL, ZSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, DCMPLX, DIMAG, SIGN
C     .. Executable Statements ..
C
      IF (N.LE.0) THEN
         TAU = ZERO
         RETURN
      END IF
C
      XNORM = DZNRM2(N-1,X,INCX)
      ALPHR = DBLE(ALPHA)
      ALPHI = DIMAG(ALPHA)
C
      IF (XNORM.EQ.ZERO .AND. ALPHI.EQ.ZERO) THEN
C
C        H  =  I
C
         TAU = ZERO
      ELSE
C
C        general case
C
         BETA = -SIGN(F08ASU(ALPHR,ALPHI,XNORM),ALPHR)
         SAFMIN = X02AMF()
         RSAFMN = ONE/SAFMIN
C
         IF (ABS(BETA).LT.SAFMIN) THEN
C
C           XNORM, BETA may be inaccurate; scale X and recompute them
C
            KNT = 0
   20       CONTINUE
            KNT = KNT + 1
            CALL ZDSCAL(N-1,RSAFMN,X,INCX)
            BETA = BETA*RSAFMN
            ALPHI = ALPHI*RSAFMN
            ALPHR = ALPHR*RSAFMN
            IF (ABS(BETA).LT.SAFMIN) GO TO 20
C
C           New BETA is at most 1, at least SAFMIN
C
            XNORM = DZNRM2(N-1,X,INCX)
            ALPHA = DCMPLX(ALPHR,ALPHI)
            BETA = -SIGN(F08ASU(ALPHR,ALPHI,XNORM),ALPHR)
            TAU = DCMPLX((BETA-ALPHR)/BETA,-ALPHI/BETA)
            ALPHA = F06CLF(DCMPLX(ONE),ALPHA-BETA,DIVFLG)
            CALL ZSCAL(N-1,ALPHA,X,INCX)
C
C           If ALPHA is subnormal, it may lose relative accuracy
C
            ALPHA = BETA
            DO 40 J = 1, KNT
               ALPHA = ALPHA*SAFMIN
   40       CONTINUE
         ELSE
            TAU = DCMPLX((BETA-ALPHR)/BETA,-ALPHI/BETA)
            ALPHA = F06CLF(DCMPLX(ONE),ALPHA-BETA,DIVFLG)
            CALL ZSCAL(N-1,ALPHA,X,INCX)
            ALPHA = BETA
         END IF
      END IF
C
      RETURN
C
C     End of F08ASV (ZLARFG)
C
      END
