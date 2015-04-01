      SUBROUTINE G05FFF(A,B,N,X,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     G05FFF generates a vector of pseudo-random Gamma deviates with
C     shape parameters A and B.
C     A = alpha
C     B = beta
C     N = number of observations
C     X = the output vector which contains Gamma deviates
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE, TWO, THREE, QUART3
      PARAMETER         (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0,
     *                  THREE=3.0D0,QUART3=0.75D0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G05FFF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  C, D, E, F, G, LASTA, P, Q, R1, R2, S, T, U, W,
     *                  XI
      INTEGER           I, IERR, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  G05CAF, X02AKF
      INTEGER           P01ABF
      EXTERNAL          G05CAF, X02AKF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         EXP, LOG, SQRT
C     .. Save statement ..
      SAVE              LASTA, T, P, S, Q, U, E, C
C     .. Data statements ..
      DATA              LASTA/-1.0D0/
C     .. Executable Statements ..
C
      IERR = 0
      NREC = 2
      IF (A.LE.ZERO .OR. B.LE.ZERO .OR. N.LE.0) THEN
         IERR = 1
         WRITE (REC,FMT=99999) A, B, N
      ELSE
C
C     The G6 algorithm
C
         IF (A.LT.ONE) THEN
            IF (LASTA.NE.A) THEN
               LASTA = A
               U = X02AKF()
               T = ONE - A
               P = T/(T+A*EXP(-T))
               Q = P*(U/T)**A
               S = ONE/A
            END IF
C
            DO 40 I = 1, N
   20          R1 = G05CAF(R1)
               IF (R1.LE.Q) THEN
                  XI = ZERO
               ELSE IF (R1.LE.P) THEN
                  XI = T*(R1/P)**S
                  W = XI
               ELSE
                  XI = T + LOG((ONE-P)/(ONE-R1))
                  W = T*LOG(XI/T)
               END IF
               R2 = G05CAF(R2)
               IF ((ONE-R2).LE.W) THEN
                  IF ((ONE/R2-ONE).LE.W) GO TO 20
                  IF ((-LOG(R2)).LE.W) GO TO 20
               END IF
               X(I) = XI*B
   40       CONTINUE
C
         ELSE IF (A.EQ.ONE) THEN
            DO 60 I = 1, N
               X(I) = -B*LOG(G05CAF(ZERO))
   60       CONTINUE
C
C     The Gbest algorithm
C
         ELSE
            IF (A.NE.LASTA) THEN
               LASTA = A
               E = A - ONE
               C = THREE*A - QUART3
            END IF
            DO 100 I = 1, N
   80          R1 = G05CAF(R1)
               G = R1 - R1*R1
               F = (R1-HALF)*(SQRT(C/G))
               XI = E + F
               IF (XI.LE.ZERO) GO TO 80
               R2 = G05CAF(R2)
               D = 64.0D0*R2*R2*G*G*G
               IF (D.GE.(ONE-TWO*F*F/XI)) THEN
                  IF (LOG(D).GE.(TWO*(E*LOG(XI/E)-F))) GO TO 80
               END IF
               X(I) = XI*B
  100       CONTINUE
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (' ** On Entry,   A, B, or N .le. 0',/'    A = ',D13.5,4X,
     *       'B = ',D13.5,4X,'N = ',I16)
      END
