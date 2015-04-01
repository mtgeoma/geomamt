      SUBROUTINE G11SAU(Y,XL,X,P,S,N2,AX,XN,Q,ALPHA,GAMMA,A,C,GPROB,PHI,
     *                  NROWXR)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     CALCULATE THE FACTOR SCORES
C
C     .. Scalar Arguments ..
      INTEGER           N2, NROWXR, Q, S
      LOGICAL           GPROB
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N2), ALPHA(N2), AX(20), C(N2), GAMMA(N2),
     *                  P(S), PHI(N2,20), XL(S), XN(20), Y(S)
      LOGICAL           X(NROWXR,N2)
C     .. Local Scalars ..
      DOUBLE PRECISION  A2, LIMIT, PROD, SUM
      INTEGER           IFAIL, J, K, L
C     .. Local Arrays ..
      DOUBLE PRECISION  LX(20)
C     .. External Functions ..
      DOUBLE PRECISION  S15ABF, X02AMF
      EXTERNAL          S15ABF, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC         EXP, LOG
C     .. Executable Statements ..
      LIMIT = LOG(X02AMF())
C
C     CALCULATE THE ELEMENTS OF PHI USING A STANDARD ROUTINE
C
      IF (GPROB) THEN
         DO 40 J = 1, N2
            DO 20 K = 1, Q
               A2 = C(J) + A(J)*XN(K)
               IFAIL = 0
               PHI(J,K) = S15ABF(A2,IFAIL)
   20       CONTINUE
   40    CONTINUE
      ELSE
         DO 120 J = 1, N2
            DO 100 K = 1, Q
               SUM = C(J) + A(J)*XN(K)
               IF (SUM.LE.0.0D0) THEN
                  IF ((SUM).LT.(LIMIT)) THEN
                     SUM = 0
                     GO TO 60
                  END IF
                  SUM = EXP(SUM)
   60             PHI(J,K) = SUM/(1.0D0+SUM)
               ELSE
                  IF ((-SUM).LT.(LIMIT)) THEN
                     PHI(J,K) = 1.0D0
                     GO TO 80
                  END IF
                  PHI(J,K) = 1.0D0/(1.0D0+EXP(-SUM))
   80          END IF
  100       CONTINUE
  120    CONTINUE
C
      END IF
C
C     CALCULATE THE LX'S
C
      DO 240 L = 1, S
C
         DO 160 K = 1, Q
            PROD = 1.0D0
            DO 140 J = 1, N2
C
               IF (X(L,J)) THEN
                  PROD = PROD*PHI(J,K)
               ELSE
                  PROD = PROD*(1.0D0-PHI(J,K))
               END IF
C
  140       CONTINUE
            LX(K) = PROD
  160    CONTINUE
C
C        CALCULATE P(L)
C
         SUM = 0.0D0
         DO 180 K = 1, Q
            SUM = SUM + LX(K)*AX(K)
  180    CONTINUE
         P(L) = SUM
C
C        CALCULATE THE FACTOR SCORES
C
         SUM = 0.0D0
         DO 200 K = 1, Q
            SUM = SUM + XN(K)*LX(K)*AX(K)
  200    CONTINUE
         Y(L) = SUM/P(L)
         IF ( .NOT. GPROB) THEN
            SUM = 0.0D0
            DO 220 J = 1, N2
               IF (X(L,J)) SUM = SUM + A(J)
  220       CONTINUE
            XL(L) = SUM
         END IF
  240 CONTINUE
C
      RETURN
      END
