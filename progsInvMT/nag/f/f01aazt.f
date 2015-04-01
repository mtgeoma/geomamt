      SUBROUTINE F01AAZ(A,LDA,N,X,LDX,P,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     APPROXIMATE INVERSE OF A REAL MATRIX
C     1ST AUGUST 1971
C
C     Originally called F01AAF
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01AAZ')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, LDX, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,N), P(N), X(LDX,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  D1, EPS
      INTEGER           I, ISAVE, J, JP
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F03AFF, DGEMV, DSCAL, DSWAP, DTRMV
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL = 1
C     Compute LU factorization of A
      EPS = X02AJF()
      CALL F03AFF(N,EPS,A,LDA,D1,I,P,IFAIL)
      IF (IFAIL.EQ.0) GO TO 20
      IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
   20 CONTINUE
C     Copy U to array X
      DO 80 J = 1, N
         DO 40 I = 1, J - 1
            X(I,J) = A(I,J)
   40    CONTINUE
         X(J,J) = 1.0D0
         DO 60 I = J + 1, N
            X(I,J) = 0.0D0
   60    CONTINUE
   80 CONTINUE
C     Compute inverse of U in array X, overwriting U
      DO 120 J = 1, N
         CALL DTRMV('U','N','U',J-1,X,LDX,X(1,J),1)
         DO 100 I = 1, J - 1
            X(I,J) = -X(I,J)
  100    CONTINUE
  120 CONTINUE
C     Compute  X * inv(L)
      DO 140 J = N, 1, -1
         IF (J.LT.N) CALL DGEMV('N',N,N-J,-1.0D0,X(1,J+1),LDX,A(J+1,J),
     *                          1,1.0D0,X(1,J),1)
         CALL DSCAL(N,1.0D0/A(J,J),X(1,J),1)
  140 CONTINUE
C     Permute columns of X
      DO 160 J = N, 1, -1
         JP = P(J) + 0.5D0
         IF (JP.NE.J) CALL DSWAP(N,X(1,J),1,X(1,JP),1)
  160 CONTINUE
      RETURN
      END
