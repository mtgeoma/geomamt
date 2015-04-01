      SUBROUTINE E04UPS(N,LDZY,NFREE,UNITQ,KX,M,FJAC,LDFJ,R,LDR,ZY,WORK)
C     MARK 14 RELEASE.  NAG COPYRIGHT 1989.
C     MARK 15A REVISED. IER-904 (APR 1991).
C     MARK 16 REVISED. IER-1099 (JUL 1993).
C
C     ******************************************************************
C     E04UPS loads the Cholesky factor of Q'HQ, where  H = J'J, the
C     approximate Hessian.
C
C     Systems Optimization Laboratory, Stanford University, California.
C     Original version written 9-May-1989.
C     This version of  E04UPS dated  9-May-1989.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           LDFJ, LDR, LDZY, M, N, NFREE
      LOGICAL           UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  FJAC(LDFJ,*), R(LDR,*), WORK(N), ZY(LDZY,*)
      INTEGER           KX(N)
C     .. Local Scalars ..
      INTEGER           I, INFO, J, L, MR, NFIXED
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, F01QCF, F01QDF, F06QHF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
      NFIXED = N - NFREE
      MR = MIN(M,N)
C
C     ------------------------------------------------------------------
C     Form the QR factorization of FJAC.
C     Note that F01QCF requires M .ge. N.
C     ------------------------------------------------------------------
      CALL F01QCF(M,MR,FJAC,LDFJ,WORK,INFO)
C
      IF (M.LT.N) THEN
         INFO = 0
         CALL F01QDF('Transpose','Separate',M,MIN(M,N),FJAC,LDFJ,WORK,
     *               N-M,FJAC(1,M+1),LDFJ,WORK(M+1),INFO)
      END IF
C
      IF (NFIXED.GT.0) THEN
         DO 40 I = 1, MR
            DO 20 L = 1, NFIXED
               J = KX(NFREE+L)
               IF (I.LE.J) THEN
                  R(I,NFREE+L) = FJAC(I,J)
               ELSE
                  R(I,NFREE+L) = ZERO
               END IF
   20       CONTINUE
   40    CONTINUE
      END IF
C
      IF (NFREE.GT.0) THEN
         DO 80 I = 1, MR
            DO 60 L = 1, NFREE
               J = KX(L)
               IF (I.LE.J) THEN
                  R(I,L) = FJAC(I,J)
               ELSE
                  R(I,L) = ZERO
               END IF
   60       CONTINUE
   80    CONTINUE
      END IF
C
      IF ((NFREE.GT.0) .AND. ( .NOT. UNITQ)) THEN
         DO 100 I = 1, MR
            CALL DGEMV('Transpose',NFREE,NFREE,ONE,ZY,LDZY,R(I,1),LDR,
     *                 ZERO,WORK,1)
            CALL DCOPY(NFREE,WORK,1,R(I,1),LDR)
  100    CONTINUE
      END IF
C
      CALL F01QCF(MR,MR,R,LDR,WORK,INFO)
C
      IF (MR.LT.N) THEN
         INFO = 0
         CALL F01QDF('Transpose','Separate',MR,MR,R,LDR,WORK,N-MR,
     *               R(1,MR+1),LDR,WORK(MR+1),INFO)
         CALL F06QHF('General',N-MR,N,ZERO,ZERO,R(MR+1,1),LDR)
      END IF
C
      RETURN
C
C     End of E04UPS.  (NLJTJ)
C
      END
