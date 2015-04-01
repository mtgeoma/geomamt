      SUBROUTINE D02RAS(A,C,DEL,Y,M,N,P,R,IR,IC,U,MTNMAX,NMAX,X,SING)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 10B REVISED. IER-401 (JAN 1983).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     -------------------------------------------------------
C
C     SOLVE  IS A COMPANION SUBROUTINE TO  DECOMP. IT SOLVES THE
C     LINEAR SYSTEM OF EQUATIONS
C
C     AA * U = X
C
C     WHERE  Y  IS A  N*M  VECTOR, AND  AA IS DESCRIBED IN
C     DECOMP.
C     SOLVE  USES THE  LU  DECOMPOSITION CONSTRUCTED IN
C     DECOMP.
C     IT PERFORMS THE NECESSARY PERMUTATIONS ON THE RIGHT
C     HAND SIDE  Y  AND ON THE ANSWER VECTOR  U.
C     FOR DETAILS ON THE STORAGE MODE OF THE SPARSE MATRIX AA
C     AND THE VARIOUS PARAMETERS SEE DECOMP.
C
C     .. Scalar Arguments ..
      INTEGER           M, MTNMAX, N, NMAX, P, R
      LOGICAL           SING
C     .. Array Arguments ..
      DOUBLE PRECISION  A(MTNMAX,M), C(MTNMAX,M), DEL(M,MTNMAX),
     *                  U(MTNMAX), X(MTNMAX), Y(MTNMAX)
      INTEGER           IC(NMAX,M), IR(NMAX,M)
C     .. Local Scalars ..
      DOUBLE PRECISION  TE
      INTEGER           I, I0, I0J, I0PK, I1, I11, I11PK, I2, I2PK, II,
     *                  IX, IXN, IXP, IXPJ, J, J1, J1MI0, JM, K, K1, M1,
     *                  MN, MNPI, N1, N1M, N2, P1, P2
C     .. Executable Statements ..
      M1 = M - 1
      N2 = N + 1
      N1 = N - 1
      P1 = P - 1
      P2 = P + 1
C
C     ROW INTERCHANGES ON THE RIGHT HAND SIDE
C
C     -------------------------------------------------------
C
      SING = .FALSE.
      IF (P.EQ.0) GO TO 40
      DO 20 I = 1, P
         X(I) = Y(I)
   20 CONTINUE
   40 MN = M*N1
      DO 60 I = 1, M
         MNPI = MN + I
         X(MNPI) = Y(MNPI)
   60 CONTINUE
      DO 100 I = 1, N1
         IX = (I-1)*M + P
         DO 80 J = 1, M
            IXP = IX + IR(I,J)
            IXPJ = IX + J
            X(IXPJ) = Y(IXP)
   80    CONTINUE
  100 CONTINUE
C
C     SOLVE  L * Y = X
C
      DO 160 I = 2, N
         I1 = M*(I-2) + 1
         I2 = I1 + P1
         I11 = I1 - 1
         IF (P.EQ.0) GO TO 180
         DO 140 J = I1, I2
            JM = J + M
            TE = X(JM)
            DO 120 K = 1, M
               I11PK = I11 + K
               TE = TE - C(J,K)*X(I11PK)
  120       CONTINUE
            X(JM) = TE
  140    CONTINUE
  160 CONTINUE
C
  180 N1M = N1*M
      IXN = N1M + P
      IF (R.EQ.0) GO TO 240
      DO 220 II = 1, R
         I1 = IXN + II
         TE = X(I1)
         DO 200 J = 1, N1M
            TE = TE - DEL(II,J)*X(J)
  200    CONTINUE
         X(I1) = TE
  220 CONTINUE
C
C     SOLVE  U * Z = Y
C
  240 DO 440 I = 1, N
         II = N2 - I
         I0 = (II-1)*M
         I1 = I0 + 2
         I2 = II*M
         IF (I.EQ.1) GO TO 300
         DO 280 J = P2, M
            I0J = I0 + J
            TE = X(I0J)
            DO 260 K = 1, M
               I2PK = I2 + K
               TE = TE - C(I0J,K)*X(I2PK)
  260       CONTINUE
            X(I0J) = TE
  280    CONTINUE
  300    DO 340 J = I1, I2
            K1 = J - I1 + 1
            TE = X(J)
            DO 320 K = 1, K1
               I0PK = I0 + K
               TE = TE - A(J,K)*X(I0PK)
  320       CONTINUE
            X(J) = TE
  340    CONTINUE
         IF (A(I2,M).NE.0.0D0) GO TO 360
         SING = .TRUE.
         RETURN
  360    X(I2) = X(I2)/A(I2,M)
         DO 420 J = 1, M1
            J1 = I2 - J
            TE = X(J1)
            K1 = M - J + 1
            DO 380 K = K1, M
               I0PK = I0 + K
               TE = TE - A(J1,K)*X(I0PK)
  380       CONTINUE
            J1MI0 = J1 - I0
            IF (A(J1,J1MI0).NE.0.0D0) GO TO 400
            SING = .TRUE.
            RETURN
  400       X(J1) = TE/A(J1,J1MI0)
  420    CONTINUE
  440 CONTINUE
C
C     INTERCHANGES IN  X
C
      DO 480 I = 1, N
         IX = (I-1)*M
         DO 460 J = 1, M
            IXP = IX + IC(I,J)
            IXPJ = IX + J
            U(IXP) = X(IXPJ)
  460    CONTINUE
  480 CONTINUE
C
      RETURN
      END
