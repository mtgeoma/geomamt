      SUBROUTINE F01APF(N,LOW,IUPP,INTGER,H,IH,V,IV)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C     MARK 5C REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     DIRTRANS
C     FORM THE MATRIX OF ACCUMULATED TRANSFORMATIONS IN THE ARRAY
C     V(N,N) FROM THE INFORMATION LEFT BY SUBROUTINE F01AKF
C     BELOW THE UPPER HESSENBERG MATRIX, H, IN THE ARRAY H(N,N)
C     AND IN THE INTEGER ARRAY INTGER(N).
C     1ST AUGUST 1971
C
C     .. Scalar Arguments ..
      INTEGER           IH, IUPP, IV, LOW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  H(IH,N), V(IV,N)
      INTEGER           INTGER(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  X
      INTEGER           I, I1, II, J, LOW1, M
C     .. Executable Statements ..
      DO 40 I = 1, N
         DO 20 J = 1, N
            V(I,J) = 0.0D0
   20    CONTINUE
         V(I,I) = 1.0D0
   40 CONTINUE
      LOW1 = LOW + 1
      IF (LOW1.GT.IUPP) RETURN
      DO 120 II = LOW1, IUPP
         I = LOW1 + IUPP - II
         I1 = I - 1
         IF (LOW1.GT.I1) GO TO 80
         DO 60 J = LOW1, I1
            V(I,J) = H(I,J-1)
   60    CONTINUE
   80    M = INTGER(I)
         IF (M.EQ.I) GO TO 120
         DO 100 J = LOW1, IUPP
            X = V(M,J)
            V(M,J) = V(I,J)
            V(I,J) = X
  100    CONTINUE
  120 CONTINUE
      RETURN
      END
