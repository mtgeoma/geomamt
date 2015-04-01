      SUBROUTINE Y90DKF(MATRIX,A,NDA,N)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ==================================================
C         *  Y90DKF :  Symmetrize a Complex Square Matrix  *
C         ==================================================
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      COMPLEX*16        CHALF
      DOUBLE PRECISION  ZERO
      PARAMETER         (CHALF=(0.5D0,0.0D0),ZERO=0.0D0)
C     .. Scalar Arguments ..
      INTEGER           N, NDA
      CHARACTER*1       MATRIX
C     .. Array Arguments ..
      COMPLEX*16        A(NDA,*)
C     .. Local Scalars ..
      COMPLEX*16        Z
      INTEGER           I, J
C     .. External Functions ..
      LOGICAL           Y90WAF
      EXTERNAL          Y90WAF
C     .. Intrinsic Functions ..
      INTRINSIC         DCMPLX, DCONJG, DBLE
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Symmetrize the matrix
C
C-----------------------------------------------------------------------
C
C     Copy upper triangular part into lower triangular part
C
      IF (Y90WAF(MATRIX,'U')) THEN
         DO 40 J = 1, N
            DO 20 I = J + 1, N
               A(I,J) = DCONJG(A(J,I))
   20       CONTINUE
   40    CONTINUE
C
C     Copy lower triangular part into upper triangular part
C
      ELSE IF (Y90WAF(MATRIX,'L')) THEN
         DO 80 J = 1, N
            DO 60 I = J + 1, N
               A(J,I) = DCONJG(A(I,J))
   60       CONTINUE
   80    CONTINUE
C
C     Average symmetric elements
C
      ELSE
         DO 120 J = 1, N
            DO 100 I = J + 1, N
               Z = CHALF*(A(I,J)+DCONJG(A(J,I)))
               A(I,J) = Z
               A(J,I) = DCONJG(Z)
  100       CONTINUE
  120    CONTINUE
      END IF
C
C     Set to zero the imaginary parts of te diagonal elements
C
      DO 140 I = 1, N
         A(I,I) = DCMPLX(DBLE(A(I,I)),ZERO)
  140 CONTINUE
C-----------------------------------------------------------------------
C
C     End of Y90DKF
C
C-----------------------------------------------------------------------
      RETURN
      END
