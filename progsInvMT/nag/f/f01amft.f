      SUBROUTINE F01AMF(N,K,L,AR,IAR,AI,IAI,INTGER)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4.5 REVISED
C     MARK 11 REVISED. VECTORISATION (JAN 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     COMHES
C     GIVEN THE COMPLEX UNSYMMETRIC MATRIX,A, STORED IN THE ARRAYS
C     AR(N,N) AND AI(N,N), THIS SUBROUTINE REDUCES THE
C     SUB-MATRIX OF ORDER L - K + 1, WHICH STARTS AT THE ELEMENTS
C     AR(K,K) AND AI(K,K) AND FINISHES AT THE ELEMENTS AR(L,L)
C     AND AI(L,L), TO HESSENBERG FORM, H, BY NON-ORTHOGONAL
C     ELEMENTARY TRANSFORMATIONS. THE MATRIX H IS OVERWRITTEN ON A
C     WITH DETAILS OF THE TRANSFORMATIONS STORED IN THE REMAINING
C     TRIANGLE UNDER H AND IN ELEMENTS K TO L OF THE ARRAY
C     INTGER(N).
C     1ST AUGUST 1971
C
C     .. Scalar Arguments ..
      INTEGER           IAI, IAR, K, L, N
C     .. Array Arguments ..
      DOUBLE PRECISION  AI(IAI,N), AR(IAR,N)
      INTEGER           INTGER(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  XI, XR, YI, YR
      INTEGER           I, J, K1, LA, M, M1, MM
C     .. External Subroutines ..
      EXTERNAL          A02ACF, F01AMY, F01AMZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      LA = L - 1
      K1 = K + 1
      IF (K1.GT.LA) RETURN
      DO 140 M = K1, N
         I = M
         XR = 0.0D0
         XI = 0.0D0
         IF (M.GT.L) GO TO 120
         DO 20 J = M, L
            IF ((ABS(AR(J,M-1))+ABS(AI(J,M-1))).LE.(ABS(XR)+ABS(XI)))
     *          GO TO 20
            XR = AR(J,M-1)
            XI = AI(J,M-1)
            I = J
   20    CONTINUE
         INTGER(M) = I
         IF (I.EQ.M) GO TO 80
C        INTERCHANGE ROWS AND COLUMNS OF ARRAYS AR AND AI.
         DO 40 J = K, N
            YR = AR(I,J)
            AR(I,J) = AR(M,J)
            AR(M,J) = YR
            YI = AI(I,J)
            AI(I,J) = AI(M,J)
            AI(M,J) = YI
   40    CONTINUE
         DO 60 J = 1, L
            YR = AR(J,I)
            AR(J,I) = AR(J,M)
            AR(J,M) = YR
            YI = AI(J,I)
            AI(J,I) = AI(J,M)
            AI(J,M) = YI
   60    CONTINUE
   80    IF (XR.EQ.0.0D0 .AND. XI.EQ.0.0D0) GO TO 120
         IF (M.GE.L) GO TO 120
         M1 = M + 1
         DO 100 I = M1, L
            YR = AR(I,M-1)
            YI = AI(I,M-1)
            CALL A02ACF(YR,YI,XR,XI,AR(I,M-1),AI(I,M-1))
  100    CONTINUE
         CALL F01AMZ(AR(1,M+1),IAR,AI(1,M+1),IAI,L,L-M,AR(M+1,M-1)
     *               ,1,AI(M+1,M-1),1,AR(1,M),AI(1,M))
  120    CALL F01AMY(AR(K+1,K),IAR,AI(K+1,K),IAI,L-K,M-K,AR(K+1,M)
     *               ,AI(K+1,M))
  140 CONTINUE
      DO 180 MM = K1, LA
         M = LA - MM + K1
         M1 = M - 2
         I = INTGER(M)
         IF (I.EQ.M) GO TO 180
         IF (K.GT.M1) GO TO 180
         DO 160 J = K, M1
            YR = AR(I,J)
            AR(I,J) = AR(M,J)
            AR(M,J) = YR
            YI = AI(I,J)
            AI(I,J) = AI(M,J)
            AI(M,J) = YI
  160    CONTINUE
  180 CONTINUE
      RETURN
      END
