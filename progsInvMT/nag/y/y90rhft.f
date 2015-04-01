      SUBROUTINE Y90RHF(DTYPE,TYPE,N,KL,NROW,A,D,DIAG,ODIAG,COND,SCALE,
     *                  DETMAN,DETEXP,DIST,SEED,IROW,WORK1,WORK2,IWORK2)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ===================================================
C         *  Y90RHF :  Generate Variable Bandwidth Matrices  *
C         ===================================================
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, HALF, TWO
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0,TWO=2.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  COND, DETMAN, SCALE
      INTEGER           DETEXP, DIST, DTYPE, IWORK2, KL, N, TYPE
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*), D(*), DIAG(*), ODIAG(*), WORK1(*),
     *                  WORK2(IWORK2,*)
      INTEGER           IROW(*), NROW(*), SEED(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  C, PI, S, T, X, X1, X2
      INTEGER           I, J, J1, J2, K1, K2, L, M, M1, TTYPE
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF, Y90TBF
      EXTERNAL          X01AAF, Y90TBF
C     .. External Subroutines ..
      EXTERNAL          Y90RAF, Y90SMF
C     .. Intrinsic Functions ..
      INTRINSIC         ANINT, COS, MAX, MIN, NINT, DBLE, SIN
C     .. Statement Functions ..
      INTEGER           LL
C     .. Statement Function definitions ..
      LL(I,J) = IROW(I+1) + J - I - 1
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Generate a generic variable bandwidth matrix
C
C-----------------------------------------------------------------------
      IF (TYPE.LE.0) THEN
C
C     Generate a banded matrix of bandwidth KL
C
         TTYPE = 0
         CALL Y90RAF('S',DTYPE,TTYPE,N,N,KL,KL,WORK2,IWORK2,D,DIAG,
     *               ODIAG,COND,SCALE,DETMAN,DETEXP,DIST,SEED)
C
C     Generate a random plane rotation
C
         K1 = NINT(Y90TBF(1,SEED)*DBLE(N-1)) + 1
   20    CONTINUE
         K2 = NINT(Y90TBF(1,SEED)*DBLE(N-1)) + 1
         IF (K2.EQ.K1) GO TO 20
         IF (K1.GT.K2) THEN
            I = K1
            K1 = K2
            K2 = I
         END IF
C
         PI = X01AAF(PI)
         T = TWO*PI*Y90TBF(1,SEED) - PI
         C = COS(T)
         S = SIN(T)
C
C     Set up variable bandwidth storage
C
         DO 40 I = 1, N
            IROW(I+1) = MAX(1,I-KL)
   40    CONTINUE
         IROW(K2+1) = IROW(K1+1)
         DO 60 I = MAX(1,K2-KL), MIN(N,K2+KL)
            IROW(I+1) = MIN(IROW(I+1),K1)
   60    CONTINUE
C
         IROW(1) = 1
         DO 80 I = 1, N
            NROW(I) = I + 1 - IROW(I+1)
            IROW(I+1) = IROW(I) + NROW(I)
   80    CONTINUE
         DO 100 I = 1, IROW(N+1) - 1
            A(I) = ZERO
  100    CONTINUE
C
C     Copy band matrix into it
C
         DO 140 I = 1, N
            DO 120 J = MAX(1,I-KL), I
               A(LL(I,J)) = WORK2(I-J+1,J)
  120       CONTINUE
  140    CONTINUE
C
C     Apply the generated random plane rotation
C
C     1. Initialise
C
         DO 160 I = 1, K2 - K1
            J1 = LL(I+K1,K1)
            IF (J1.GE.IROW(I+K1)) THEN
               WORK1(I) = A(LL(I+K1,K1))
            ELSE
               WORK1(I) = ZERO
            END IF
  160    CONTINUE
C
C     2. Right rotation (columns)
C
         WORK1(K2-K1) = C*WORK1(K2-K1) - S*A(LL(K1,K1))
C
         DO 180 I = K1, K2 - 1
            J1 = LL(I,K1)
            IF (J1.GE.IROW(I)) A(J1) = C*A(J1) + S*A(LL(K2,I))
  180    CONTINUE
C
         DO 200 I = K2, MIN(N,K2+KL)
            J1 = LL(I,K1)
            J2 = LL(I,K2)
            X1 = A(J1)
            X2 = A(J2)
            A(J1) = C*X1 + S*X2
            A(J2) = C*X2 - S*X1
  200    CONTINUE
C
C     3. Left rotation (rows)
C
         DO 220 I = MAX(1,K1-KL), K1
            J1 = LL(K1,I)
            J2 = LL(K2,I)
            X1 = A(J1)
            X2 = A(J2)
            A(J1) = C*X1 + S*X2
            A(J2) = C*X2 - S*X1
  220    CONTINUE
C
         DO 240 I = 1, K2 - K1
            A(LL(K2,I+K1)) = C*A(LL(K2,I+K1)) - S*WORK1(I)
  240    CONTINUE
C-----------------------------------------------------------------------
C
C     Generate a (non-generic) triangular variable-bandwidth matrix
C
C-----------------------------------------------------------------------
      ELSE
         IROW(1) = 1
         L = 0
         DO 280 I = 1, N
            M1 = NINT(Y90TBF(1,SEED)*DBLE(KL+1)+0.5D0)
            M = MIN(M1,I)
            NROW(I) = M
            IROW(I+1) = IROW(I) + M
            DO 260 J = 1, M - 1
               IF (DIST.LE.1) THEN
                  X = ODIAG(1) + (ODIAG(2)-ODIAG(1))*Y90TBF(DIST,SEED)
               ELSE IF (DIST.EQ.2) THEN
                  X = HALF*((ODIAG(1)+ODIAG(2))+(ODIAG(2)-ODIAG(1))
     *                *Y90TBF(DIST,SEED))
               ELSE
                  X = DIAG(1) + DIAG(2)*Y90TBF(DIST,SEED)
               END IF
               IF (TYPE.EQ.1) X = ANINT(X)
               L = L + 1
               A(L) = X
  260       CONTINUE
            IF (DIST.LE.1) THEN
               X = DIAG(1) + (DIAG(2)-DIAG(1))*Y90TBF(DIST,SEED)
            ELSE IF (DIST.EQ.2) THEN
               X = HALF*((DIAG(1)+DIAG(2))+(DIAG(2)-DIAG(1))
     *             *Y90TBF(DIST,SEED))
            ELSE
               X = DIAG(1) + DIAG(2)*Y90TBF(DIST,SEED)
            END IF
            IF (TYPE.EQ.1) X = ANINT(X)
            L = L + 1
            A(L) = X
  280    CONTINUE
C
      END IF
C-----------------------------------------------------------------------
C
C     Calculate the determinant
C
C-----------------------------------------------------------------------
      DETMAN = ONE
      DETEXP = 0
      DO 300 I = 1, N
         DETMAN = DETMAN*A(IROW(I+1)-1)
         CALL Y90SMF(DETMAN,DETEXP,4)
  300 CONTINUE
C-----------------------------------------------------------------------
C
C     End of Y90RHF
C
C-----------------------------------------------------------------------
      RETURN
      END
