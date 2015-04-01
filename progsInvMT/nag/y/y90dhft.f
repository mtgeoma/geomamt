      SUBROUTINE Y90DHF(MATRIX,M,N,A,NDA,B,NDB)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         =====================================
C         *  Y90DHF :  Copy a Complex Matrix  *
C         =====================================
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           M, N, NDA, NDB
      CHARACTER*1       MATRIX
C     .. Array Arguments ..
      COMPLEX*16        A(NDA,*), B(NDB,*)
C     .. Local Scalars ..
      INTEGER           I, J
C     .. External Functions ..
      LOGICAL           Y90WAF
      EXTERNAL          Y90WAF
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Copy the Matrix Transpose
C
C-----------------------------------------------------------------------
      IF (Y90WAF(MATRIX,'T')) THEN
         DO 40 I = 1, N
            DO 20 J = 1, M
               B(I,J) = A(J,I)
   20       CONTINUE
   40    CONTINUE
C-----------------------------------------------------------------------
C
C     Copy the Matrix Hermitian Conjugate
C
C-----------------------------------------------------------------------
      ELSE IF (Y90WAF(MATRIX,'C')) THEN
         DO 80 I = 1, N
            DO 60 J = 1, M
               B(I,J) = DCONJG(A(J,I))
   60       CONTINUE
   80    CONTINUE
C-----------------------------------------------------------------------
C
C     Copy the Matrix
C
C-----------------------------------------------------------------------
      ELSE
         DO 120 I = 1, N
            DO 100 J = 1, M
               B(J,I) = A(J,I)
  100       CONTINUE
  120    CONTINUE
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90DHF
C
C-----------------------------------------------------------------------
      RETURN
      END
