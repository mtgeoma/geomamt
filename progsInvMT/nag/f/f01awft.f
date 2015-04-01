      SUBROUTINE F01AWF(N,LOW,LHI,M,D,ZR,IZR,ZI,IZI)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C     MARK 5C REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     COMBALBAK
C     BACKWARD TRANSFORMATION OF A SET OF RIGHT-HAND EIGENVECTORS
C     OF A BALANCED MATRIX INTO THE EIGENVECTORS OF THE ORIGINAL
C     MATRIX FROM WHICH THE BALANCED MATRIX WAS DERIVED BY A CALL
C     OF SUBROUTINE F01AVF.
C     DECEMBER 1ST.,1971
C
C     .. Scalar Arguments ..
      INTEGER           IZI, IZR, LHI, LOW, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(N), ZI(IZI,M), ZR(IZR,M)
C     .. Local Scalars ..
      DOUBLE PRECISION  S
      INTEGER           I, II, J, K, LHI1, LOW1
C     .. Executable Statements ..
      IF (LOW.GT.LHI) GO TO 60
      DO 40 I = LOW, LHI
         S = D(I)
C        LEFT-HAND EIGENVECTORS ARE BACK-TRANSFORMED IF THE FOREGOING
C        STATEMENT IS REPLACED BY S=1/D(I)
         DO 20 J = 1, M
            ZR(I,J) = ZR(I,J)*S
            ZI(I,J) = ZI(I,J)*S
   20    CONTINUE
   40 CONTINUE
   60 I = LOW
      LOW1 = LOW - 1
      IF (LOW1.LT.1) GO TO 120
      DO 100 II = 1, LOW1
         I = I - 1
         K = D(I)
         IF (K.EQ.I) GO TO 100
         DO 80 J = 1, M
            S = ZR(I,J)
            ZR(I,J) = ZR(K,J)
            ZR(K,J) = S
            S = ZI(I,J)
            ZI(I,J) = ZI(K,J)
            ZI(K,J) = S
   80    CONTINUE
  100 CONTINUE
  120 LHI1 = LHI + 1
      IF (LHI1.GT.N) RETURN
      DO 160 I = LHI1, N
         K = D(I)
         IF (K.EQ.I) GO TO 160
         DO 140 J = 1, M
            S = ZR(I,J)
            ZR(I,J) = ZR(K,J)
            ZR(K,J) = S
            S = ZI(I,J)
            ZI(I,J) = ZI(K,J)
            ZI(K,J) = S
  140    CONTINUE
  160 CONTINUE
      RETURN
      END
