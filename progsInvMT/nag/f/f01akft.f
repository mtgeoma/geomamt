      SUBROUTINE F01AKF(N,K,L,A,IA,INTGER)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C     MARK 8 REVISED. IER-248 (JUN 1980).
C     MARK 11 REVISED. VECTORISATION (JAN 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12 REVISED. EXTENDED BLAS (JUNE 1986)
C
C     DIRHES
C     AUGUST 1ST, 1971 .
C     GIVEN THE UNSYMMETRIC MATRIX, A, STORED IN THE ARRAY A(N,N),
C     THIS SUBROUTINE REDUCES THE SUB-MATRIX OF ORDER L - K + 1,
C     WHICH STARTS AT THE ELEMENT A(K,K) AND FINISHES AT THE
C     ELEMENT A(L,L), TO HESSENBERG FORM, H, BY THE DIRECT
C     METHOD(AN = NH). THE MATRIX H IS OVERWRITTEN ON A WITH
C     DETAILS OF THE TRANSFORMATIONS (N) STORED IN THE REMAINING
C     TRIANGLE UNDER H AND IN ELEMENTS K TO L OF THE ARRAY
C     INTGER(N).
C     1ST AUGUST 1971
C
C     .. Scalar Arguments ..
      INTEGER           IA, K, L, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N)
      INTEGER           INTGER(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  X, Y
      INTEGER           I, J, K1, M
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DTRSV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      K1 = K + 1
      IF (K1.GT.L) RETURN
      DO 140 J = K1, N
         M = J
         X = 0.0D0
         IF (J.GT.L) GO TO 120
         DO 20 I = J, L
            IF (ABS(A(I,J-1)).LE.ABS(X)) GO TO 20
            X = A(I,J-1)
            M = I
   20    CONTINUE
         INTGER(J) = M
         IF (M.EQ.J) GO TO 80
C        INTERCHANGE ROWS AND COLUMNS OF A.
         DO 40 I = K, N
            Y = A(M,I)
            A(M,I) = A(J,I)
            A(J,I) = Y
   40    CONTINUE
         DO 60 I = 1, L
            Y = A(I,M)
            A(I,M) = A(I,J)
            A(I,J) = Y
   60    CONTINUE
   80    IF (X.NE.0.0D0 .AND. J.LT.L) THEN
            DO 100 I = J + 1, L
               A(I,J-1) = A(I,J-1)/X
  100       CONTINUE
            CALL DGEMV('N',L,L-J,1.0D0,A(1,J+1),IA,A(J+1,J-1),1,1.0D0,
     *                 A(1,J),1)
         END IF
  120    CALL DTRSV('L','N','U',J-K,A(K+1,K),IA,A(K+1,J),1)
         IF (J.LT.L) CALL DGEMV('N',L-J,J-K,-1.0D0,A(J+1,K),IA,A(K+1,J),
     *                          1,1.0D0,A(J+1,J),1)
  140 CONTINUE
      RETURN
      END
