      SUBROUTINE G05FEF(A,B,N,X,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     G05FEF generates a vector of pseudo-random Beta deviates with
C     shape parameters A and B.
C     A = alpha
C     B = beta
C     N = number of observations
C     X = the output vector which contains Beta deviates
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, QUART, HALF, ONE, TWO, LN4, OPLN5
      PARAMETER         (ZERO=0.0D0,QUART=0.25D0,HALF=0.5D0,ONE=1.0D0,
     *                  TWO=2.0D0,LN4=1.38629436492919921875D0,
     *                  OPLN5=2.60943794250488281250D0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G05FEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, B1, BIG, D, E, G, H, K, K1, K2, K3, L,
     *                  LASTA, LASTB, M, P, R, R1, R2, RR1, RR2, S,
     *                  SMALL, T, V, W, Y, Z, ZMONE
      INTEGER           I, IERR, NREC
      LOGICAL           SWITCH
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  G05CAF, X02AKF, X02ALF
      INTEGER           P01ABF
      EXTERNAL          G05CAF, X02AKF, X02ALF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         EXP, LOG, SQRT
C     .. Save statement ..
      SAVE              LASTA, LASTB, A1, B1, G, K1, K2, K3, L, M
C     .. Data statements ..
      DATA              LASTA, LASTB/-1.0D0, -1.0D0/
C     .. Executable Statements ..
C
      IERR = 0
      NREC = 2
      IF (A.LE.ZERO .OR. B.LE.ZERO .OR. N.LE.0) THEN
         IERR = 1
         WRITE (REC,FMT=99999) A, B, N
      ELSE
C
C        The Rejection method ----- Johnk
C        Two random numbers and two exponentiation of type y exp x are
C        required.
C
         IF (A.LT.HALF .AND. B.LT.HALF) THEN
            DO 40 I = 1, N
   20          R1 = G05CAF(R1)
               R2 = G05CAF(R2)
               RR1 = R1**(1/A)
               RR2 = R2**(1/B)
               K = (RR1+RR2)
               IF (K.GT.ONE .OR. K.LE.ZERO) GO TO 20
               X(I) = RR1/K
   40       CONTINUE
         ELSE
            D = X02AKF()
            E = X02ALF()
            IF (B.LT.A) THEN
               SMALL = B
               BIG = A
               SWITCH = .TRUE.
            ELSE
               SMALL = A
               BIG = B
               SWITCH = .FALSE.
            END IF
C
C        The Rejection method ----- BB
C        The LOG Logistic envelope
C
            IF (SMALL.GT.ONE) THEN
               IF (LASTA.NE.A .OR. LASTB.NE.B) THEN
                  LASTA = A
                  LASTB = B
                  A1 = SMALL + BIG
                  B1 = SQRT((A1-TWO)/(TWO*SMALL*BIG-A1))
                  G = SMALL + 1/B1
               END IF
               DO 80 I = 1, N
   60             R1 = G05CAF(R1)
                  R2 = G05CAF(R2)
                  V = B1*LOG(R1/(ONE-R1))
                  W = SMALL*EXP(V)
                  Z = R1*R1*R2
                  R = G*V - LN4
                  S = SMALL + R - W
C
                  IF (S+OPLN5.LT.5.0D0*Z) THEN
                     T = LOG(Z)
                     IF (S.LT.T) THEN
                        IF (R+A1*LOG(A1/(BIG+W)).LT.T) GO TO 60
                     END IF
                  END IF
                  IF (SWITCH) THEN
                     X(I) = BIG/(BIG+W)
                  ELSE
                     X(I) = W/(BIG+W)
                  END IF
   80          CONTINUE
C
C        The Switching method ----- Atkinson
C        Each trial requires the generation of two random numbers and
C        the equivalent of four LOG evaluations.
C
            ELSE IF (BIG.GT.ONE .AND. SMALL.LT.ONE) THEN
               IF (LASTA.NE.A .OR. LASTB.NE.B) THEN
                  LASTA = A
                  LASTB = B
                  L = (ONE-SMALL)/(BIG+ONE-SMALL)
                  M = BIG*L/((BIG*L+SMALL*(ONE-L)**BIG))
               END IF
               DO 120 I = 1, N
  100             R1 = G05CAF(R1)
                  R2 = G05CAF(R2)
                  IF (R1.LE.M) THEN
                     Y = L*((R1/M)**(ONE/SMALL))
                     IF (LOG(R2).GT.((BIG-ONE)*LOG(ONE-Y))) GO TO 100
                  ELSE
                     Y = ONE - (ONE-L)*(((ONE-R1)/(ONE-M))**(ONE/BIG))
                     IF (LOG(R2).GT.((SMALL-ONE)*LOG(Y/L))) GO TO 100
                  END IF
                  IF (SWITCH) THEN
                     X(I) = ONE - Y
                  ELSE
                     X(I) = Y
                  END IF
  120          CONTINUE
            ELSE
C
C        The Rejection method ----- BC
C
               IF (LASTA.NE.A .OR. LASTB.NE.B) THEN
                  LASTA = A
                  LASTB = B
                  A1 = BIG + SMALL
                  B1 = 1/SMALL
                  G = 1 + BIG - SMALL
                  K1 = G*(0.0138889D0+0.0416667D0*SMALL)
     *                 /(BIG*B1-0.777778D0)
                  K2 = QUART + (HALF+QUART/G)*SMALL
                  K3 = ONE/(ONE+(BIG/(SMALL*E))**SMALL)
               END IF
C
               DO 200 I = 1, N
  140             R1 = G05CAF(R1)
                  R2 = G05CAF(R2)
                  P = R1*R2
                  Z = R1*P
                  IF (R1.LT.HALF) THEN
                     IF (R1.LT.D .OR. Z.LE.ZERO) GO TO 140
                     IF ((QUART*R2+Z-P).GE.K1) GO TO 140
                  ELSE
                     IF (Z.GE.K2) GO TO 140
                     IF (Z.LE.QUART) THEN
                        V = B1*LOG(R1/(ONE-R1))
                        W = BIG*EXP(V)
                        GO TO 180
                     ELSE
                        GO TO 160
                     END IF
                  END IF
                  IF (R1.GE.K3) THEN
                     IF (4.0D0*Z.GT.(ONE+SMALL/BIG)**A1) GO TO 140
                     IF (SWITCH) THEN
                        X(I) = ONE
                     ELSE
                        X(I) = ZERO
                     END IF
                     GO TO 200
                  END IF
  160             V = B1*LOG(R1/(ONE-R1))
                  W = BIG*EXP(V)
                  H = A1*V + A1*LOG((ONE+BIG/SMALL)/(ONE+W)) - LN4
                  ZMONE = Z - ONE
                  IF (ZMONE.GT.H) THEN
                     IF (ZMONE.GT.Z*H) GO TO 140
                     IF (LOG(Z).GT.H) GO TO 140
                  END IF
C
  180             IF (SWITCH) THEN
                     X(I) = W/(SMALL+W)
                  ELSE
                     X(I) = SMALL/(SMALL+W)
                  END IF
  200          CONTINUE
            END IF
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (' ** On Entry,  A, B, or N .le. 0',/'     A = ',D13.5,4X,
     *       'B = ',D13.5,4X,'N = ',I16)
      END
