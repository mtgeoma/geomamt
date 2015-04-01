      SUBROUTINE F01BLZ(N,A,IA,P,IFAIL)
C     NAG COPYRIGHT 1976. MARK  5 RELEASE.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     PDK, DNAC, NPL, TEDDINGTON. NOVEMBER 1975.
C     NPL DNAC LIBRARY SUBROUTINE UDUT.
C
C     A CHOLESKY TYPE DECOMPOSITION IS GIVEN FOR A SYMMETRIC
C     POSITIVE DEFINITE MATRIX A OF ORDER N, SUCH THAT A=UDL
C     WHERE L IS THE TRANSPOSE OF U, WHICH IS A UNIT UPPER-
C     TRIANGULAR MATRIX,  AND D IS A DIAGONAL MATRIX. THE
C     LOWER TRIANGLE OF A IS SUPPLIED IN THE ARRAY A(IA,N)
C     AND U, OMITTING THE UNIT DIAGONAL, IS FORMED IN THE
C     REMAINING STRICT UPPER TRIANGLE. RECIPROCALS OF THE
C     DIAGONAL ELEMENTS OF D ARE STORED IN THE ARRAY P(N). A
C     FAILURE EXIT WILL OCCUR IF A, MODIFIED BY THE ROUNDING
C     ERRORS, IS NOT POSITIVE DEFINITE.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01BLZ')
C     .. Scalar Arguments ..
      INTEGER           IA, IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), P(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ONE, X, Y, Z, ZERO
      INTEGER           I, II, ISAVE, J, JJ, JP1, K
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Data statements ..
      DATA              ZERO/0.0D0/, ONE/1.0D0/
C     .. Executable Statements ..
      ISAVE = IFAIL
      DO 160 II = 1, N
         I = N + 1 - II
         DO 140 JJ = I, N
            J = N + I - JJ
            X = A(J,I)
            JP1 = J + 1
            IF (I.NE.J) GO TO 80
            IF (J.EQ.N) GO TO 40
            DO 20 K = JP1, N
               Y = A(I,K)
               Z = Y*P(K)
               A(I,K) = Z
               X = X - Y*Z
   20       CONTINUE
   40       IF (X.NE.ZERO) GO TO 60
            IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
            RETURN
   60       P(I) = ONE/X
            GO TO 140
   80       IF (J.EQ.N) GO TO 120
            DO 100 K = JP1, N
               X = X - A(I,K)*A(J,K)
  100       CONTINUE
  120       A(I,J) = X
  140    CONTINUE
  160 CONTINUE
      IFAIL = 0
      RETURN
      END
