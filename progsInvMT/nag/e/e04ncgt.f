      SUBROUTINE E04NCG(TASK,UNITQ,NFREE,N,NRANK,LDQ,LDR,KX,R,Q,V,W)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     ==================================================================
C     E04NCG forms the Cholesky factor of the Hessian.
C
C     Systems Optimization Laboratory, Stanford University.
C     Mathematics Department,          UC San Diego.
C     Original version written by PEG, 07-Jul-94.
C     This version of E04NCG dated 08-Jul-94.
C     ==================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           LDQ, LDR, N, NFREE, NRANK
      LOGICAL           UNITQ
      CHARACTER         TASK
C     .. Array Arguments ..
      DOUBLE PRECISION  Q(LDQ,*), R(LDR,*), V(N), W(N)
      INTEGER           KX(N)
C     .. Local Scalars ..
      INTEGER           I, INFO, J, M, NZ
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, E04NBW, F01QCF, F01QDF, F06FBF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
      IF (TASK.EQ.'H') THEN
C        ---------------------------------------------------------------
C        Form the triangular factor of the Hessian.
C        ---------------------------------------------------------------
C        First,  form the square matrix  R  such that  P'HP = R'R,
C        where P is the permutation  kx.
C        Compute the  QR  factorization of  R.
C
         DO 20 J = 1, N
            IF (J.GT.1) CALL F06FBF(J-1,ZERO,V,1)
            CALL DCOPY(N-J+1,R(J,J),LDR,V(J),1)
            CALL E04NBW(3,N,NZ,NFREE,LDQ,UNITQ,KX,V,Q,W)
            CALL DCOPY(N,V,1,R(J,1),LDR)
   20    CONTINUE
C
         CALL F01QCF(N,N,R,LDR,W,INFO)
C
      ELSE IF (TASK.EQ.'P') THEN
C        ---------------------------------------------------------------
C        Form the factor of the permuted Hessian  P'HP.
C        ---------------------------------------------------------------
         IF (UNITQ) THEN
C           Relax, nothing needs to be done.
         ELSE
C
            M = MIN(NRANK,NFREE)
            DO 40 I = 1, M
C
C              Set  v' = (ith row of R)*Q'.
C
               CALL DGEMV('No Transpose',NFREE,NFREE-I+1,ONE,Q(1,I),LDQ,
     *                    R(I,I),LDR,ZERO,V,1)
               CALL DCOPY(NFREE,V,1,R(I,1),LDR)
   40       CONTINUE
C
            CALL F01QCF(M,M,R,LDR,W,INFO)
C
            IF (M.LT.N) THEN
               INFO = 0
               CALL F01QDF('Transpose','Separate',M,M,R,LDR,W,N-M,
     *                     R(1,M+1),LDR,W(M+1),INFO)
            END IF
         END IF
      END IF
C
C     For safety, zero out the lower-triangular part of R.
C
      DO 60 J = 1, N - 1
         CALL F06FBF(N-J,ZERO,R(J+1,J),1)
   60 CONTINUE
C
C     End of E04NCG. (LSFRMH)
C
      RETURN
      END
