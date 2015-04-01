      DOUBLE PRECISION FUNCTION G01JDZ(M,A,C)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Computes probability using Pan's procedure
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 C
      INTEGER                          M
C     .. Array Arguments ..
      DOUBLE PRECISION                 A(M)
C     .. Local Scalars ..
      DOUBLE PRECISION                 NUM, PI2, PIN, PROD, SGN, SUM,
     *                                 SUM1, U, UFLOW, V, Y
      INTEGER                          D, H, I, IS, J1, J2, J3, J4, K,
     *                                 L, N, N2, NU
C     .. External Functions ..
      DOUBLE PRECISION                 X01AAF, X02AMF
      EXTERNAL                         X01AAF, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, COS, EXP, LOG, MOD, DBLE,
     *                                 SQRT
C     .. Executable Statements ..
C
      UFLOW = LOG(X02AMF())
      PI2 = X01AAF(0.0D0)/2.0D0
      IF (M.LT.5) THEN
         N = 200
      ELSE
         N = 24
      END IF
C
      K = 1
      H = K
      I = M
      DO 20 NU = 1, M
         IF (A(NU).GE.0.0D0) GO TO 40
   20 CONTINUE
C
C     START
C
   40 NU = NU - 1
      H = M - NU
      IF (C.EQ.0.0D0) THEN
         Y = DBLE(H) - DBLE(NU)
      ELSE
         Y = C*(A(1)-A(M))
      END IF
      IF (Y.GE.0.0D0) THEN
         D = 2
         H = NU
         K = -K
         J1 = 0
         J2 = 2
         J3 = 3
         J4 = 1
      ELSE
         D = -2
         NU = NU + 1
         J1 = M - 2
         J2 = M - 1
         J3 = M + 1
         J4 = M
      END IF
      PIN = PI2/DBLE(N)
      SUM = 0.5D0*(DBLE(K)+1.0D0)
      SGN = DBLE(K)/DBLE(N)
      N2 = N + N - 1
C
      IF (MOD(H,2).EQ.0) THEN
         IS = 0
      ELSE
         IS = 1
      END IF
      DO 140 H = IS, 0, -1
         DO 120 L = J2, NU, D
            SUM1 = A(J4)
            IF (L.EQ.0) THEN
               PROD = 0.0D0
            ELSE
               PROD = A(L)
            END IF
            U = 0.5D0*(SUM1+PROD)
            V = 0.5D0*(SUM1-PROD)
            SUM1 = 0.0D0
            DO 100 I = 1, N2, 2
               Y = U - V*COS(DBLE(I)*PIN)
               NUM = Y
               IF (-C/NUM.GE.UFLOW) THEN
                  PROD = EXP(-C/NUM)
                  DO 60 K = 1, J1
                     PROD = PROD*(NUM/(Y-A(K)))
   60             CONTINUE
                  DO 80 K = J3, M
                     PROD = PROD*(NUM/(Y-A(K)))
   80             CONTINUE
                  SUM1 = SUM1 + SQRT(ABS(PROD))
               END IF
  100       CONTINUE
            SGN = -SGN
            SUM = SUM + SGN*SUM1
            J1 = J1 + D
            J3 = J3 + D
            J4 = J4 + D
  120    CONTINUE
         IF (D.EQ.2) THEN
            J3 = J3 - 1
         ELSE
            J1 = J1 + 1
         END IF
         NU = 0
         J2 = NU
  140 CONTINUE
      G01JDZ = SUM
      RETURN
      END
