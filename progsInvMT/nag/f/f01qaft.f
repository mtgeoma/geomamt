      SUBROUTINE F01QAF(M,N,A,NRA,C,NRC,Z,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 14 REVISED. IER-730 (DEC 1989).
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (HOUSQU)
C
C     F01QAF RETURNS, IN C AND Z, THE HOUSEHOLDER QU
C     FACTORIZATION OF THE M*N (M.GE.N) MATRIX A. THAT IS A
C     FACTORIZED AS
C
C     A = Q(U) , M.GT.N,   A = QU , M.EQ.N,
C          (0)
C
C     WHERE Q IS AN M*M ORTHOGONAL MATRIX AND U IS A SQUARE UPPER
C     TRIANGULAR MATRIX.
C
C     INPUT PARAMETERS.
C
C     M     - THE NUMBER OF ROWS OF A. M MUST BE AT LEAST N.
C
C     N     - THE NUMBER OF COLUMNS OF A. N MUST BE AT LEAST UNITY.
C
C     A     - THE M*N MATRIX TO BE FACTORIZED.
C
C     NRA   - ROW DIMENSION OF A AS DECLARED IN THE CALLING PROGRAM.
C             NRA MUST BE AT LEAST M.
C
C     NRC   - ROW DIMENSION OF C AS DECLARED IN THE CALLING PROGRAM.
C             NRC MUST BE AT LEAST M.
C
C     IFAIL - THE USUAL FAILURE PARAMETER. IF IN DOUBT SET
C             IFAIL TO ZERO BEFORE CALLING THIS ROUTINE.
C
C     OUTPUT PARAMETERS.
C
C     C     - AN M*N MATRIX CONTAINING DETAILS OF THE QU
C             FACTORIZATION. U IS RETURNED IN THE UPPER
C             TRIANGULAR PART OF C. THE SUB-DIAGONAL
C             ELEMENTS OF THE J(TH) COLUMN OF C CONTAIN
C             ELEMENTS Y(J+1),Y(J+2),...,Y(M) OF THE VECTOR
C             Y SUCH THAT (I-(1/Y(J))*Y*(Y**T)) IS THE
C             J(TH)  HOUSEHOLDER TRANSFORMATION MATRIX.
C             THE ROUTINE MAY BE CALLED WITH C=A.
C
C     Z     - N ELEMENT VECTOR.
C             Z(J) CONTAINS THE ELEMENT Y(J) OF THE J(TH)
C             HOUSEHOLDER TRANSFORMATION MATRIX.
C
C     IFAIL - ON NORMAL RETURN IFAIL WILL BE ZERO.
C             IF AN INPUT PARAMETER IS INCORRECTLY SUPPLIED
C             THEN IFAIL IS SET TO UNITY. NO OTHER FAILURE
C             IS POSSIBLE.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01QAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M, N, NRA, NRC
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NRA,*), C(NRC,*), Z(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIG, SMALL, TINY
      INTEGER           I, IERR, J, K, KLAST, KP1, NR
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      INTEGER           P01ABF
      EXTERNAL          X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01QAY, F01QAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MIN, SQRT
C     .. Executable Statements ..
      IERR = IFAIL
      IF (IERR.EQ.0) IFAIL = 1
C
      IF (NRA.LT.M .OR. NRC.LT.M .OR. M.LT.N .OR. N.LT.1) GO TO 100
C
      SMALL = X02AMF()
      TINY = SQRT(SMALL)
      BIG = 1.0D0/SMALL
C
      DO 40 J = 1, N
         DO 20 I = 1, M
            C(I,J) = A(I,J)
   20    CONTINUE
   40 CONTINUE
C
      IFAIL = 0
      IF (M.EQ.1) RETURN
C
      KLAST = MIN(M-1,N)
      NR = M + 1
      DO 80 K = 1, KLAST
         NR = NR - 1
C
         CALL F01QAY(NR,C(K,K),.FALSE.,Z(K),SMALL,TINY,BIG)
C
         IF (K.EQ.N) RETURN
C
         KP1 = K + 1
         DO 60 J = KP1, N
C
            CALL F01QAZ(NR,C(K,K),Z(K),C(K,J))
C
   60    CONTINUE
   80 CONTINUE
C
      RETURN
C
  100 IFAIL = P01ABF(IERR,IFAIL,SRNAME,0,P01REC)
      RETURN
      END
