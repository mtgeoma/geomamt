      SUBROUTINE F11GBP(RESID,USEXC,ITERM,N,ZETAB,BPIGAM,BETA1,BETA2,
     *                  RHO1,RHO2,X,V1,V2,W1,U,V)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-----------------------------------------------------------------------
C
C     F11GBP - Auxiliary to compute the residual vector (Lanczos
C              Method (SYMMLQ))
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BETA1, BETA2, BPIGAM, RHO1, RHO2, ZETAB
      INTEGER           ITERM, N, RESID
      LOGICAL           USEXC
C     .. Array Arguments ..
      DOUBLE PRECISION  U(N), V(N), V1(N), V2(N), W1(N), X(N)
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DSCAL
C     .. Executable Statements ..
      IF (RESID.LE.0) THEN
         CALL DCOPY(N,X,1,U,1)
         CALL DCOPY(N,V2,1,V,1)
         IF (USEXC) THEN
            CALL DAXPY(N,ZETAB,W1,1,U,1)
            CALL DSCAL(N,BPIGAM,V,1)
         ELSE
            CALL DSCAL(N,(RHO2/BETA2),V,1)
            CALL DAXPY(N,(RHO1/BETA1),V1,1,V,1)
         END IF
      ELSE IF (RESID.EQ.1) THEN
         IF (USEXC) THEN
            CALL DCOPY(N,V2,1,V,1)
            CALL DSCAL(N,BPIGAM,V,1)
         ELSE
            CALL DSCAL(N,(RHO2/BETA2),V,1)
            CALL DCOPY(N,X,1,U,1)
         END IF
      END IF
C
C     End of subroutine F11GBP
C
      RETURN
      END
