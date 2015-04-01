      SUBROUTINE F01ANF(KL,L,IR,AR,IAR,AI,IAI,INTGER,ZR,IZR,ZI,IZI,N)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4.5 REVISED
C     MARK 11 REVISED. VECTORISATION (JAN 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     COMBAK
C     GIVEN IR EIGENVECTORS, STORED AS COLUMNS IN THE ARRAYS
C     ZR(N,IR)
C     AND ZI(N,IR), OF THE COMPLEX HESSENBERG MATRIX, H,
C     WHICH WAS FORMED BY F01AMF, AND DETAILS OF THE
C     TRANSFORMATIONS AS LEFT BY F01AMF BELOW H AND IN ELEMENTS
C     K TO L OF THE ARRAY INTGER(N), THIS SUBROUTINE FORMS THE
C     EIGENVECTORS OF THE MATRIX A AND OVERWRITES THEM ON THE
C     GIVEN VECTORS.
C     1ST AUGUST 1971
C
C     .. Scalar Arguments ..
      INTEGER           IAI, IAR, IR, IZI, IZR, KL, L, N
C     .. Array Arguments ..
      DOUBLE PRECISION  AI(IAI,N), AR(IAR,N), ZI(IZI,IR), ZR(IZR,IR)
      INTEGER           INTGER(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  XI, XR
      INTEGER           I, J, K, L1, M, M1, MM
C     .. External Subroutines ..
      EXTERNAL          F01ANZ
C     .. Executable Statements ..
      K = KL + 1
      L1 = L - 1
      IF (K.GT.L1) RETURN
      IF (IR.LE.0) RETURN
      DO 40 M = K, L1
         M1 = M - 2
         I = INTGER(M)
         IF (I.EQ.M) GO TO 40
         IF (KL.GT.M1) GO TO 40
         DO 20 J = KL, M1
            XR = AR(I,J)
            AR(I,J) = AR(M,J)
            AR(M,J) = XR
            XI = AI(I,J)
            AI(I,J) = AI(M,J)
            AI(M,J) = XI
   20    CONTINUE
   40 CONTINUE
      DO 60 J = 1, IR
         CALL F01ANZ(AR(K,K-1),IAR,AI(K,K-1),IAI,L-K+1,ZR(K,J),ZI(K,J))
   60 CONTINUE
      DO 120 MM = K, L1
         M = L1 - MM + K
         I = INTGER(M)
         IF (I.EQ.M) GO TO 120
         DO 80 J = 1, IR
            XR = ZR(I,J)
            ZR(I,J) = ZR(M,J)
            ZR(M,J) = XR
            XI = ZI(I,J)
            ZI(I,J) = ZI(M,J)
            ZI(M,J) = XI
   80    CONTINUE
         M1 = M - 2
         IF (KL.GT.M1) GO TO 120
         DO 100 J = KL, M1
            XR = AR(I,J)
            AR(I,J) = AR(M,J)
            AR(M,J) = XR
            XI = AI(I,J)
            AI(I,J) = AI(M,J)
            AI(M,J) = XI
  100    CONTINUE
  120 CONTINUE
      RETURN
      END
