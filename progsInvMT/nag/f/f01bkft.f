      SUBROUTINE F01BKF(M,N,T,A,IA,D,AIJMAX,INC,IRANK,IFAIL)
C     NAG COPYRIGHT 1976. MARK  5 RELEASE.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     PDK, DNAC, NPL, TEDDINGTON. NOVEMBER 1975.
C     NPL DNAC LIBRARY SUBROUTINE QUFAC.
C
C     DETERMINE THE FACTORS OF THE HOUSEHOLDER DECOMPOSITION
C     WITH PIVOTING OF AN M*N MATRIX A AND THE RANK OF THE
C     MATRIX, WHERE M.GE.N. THE MATRIX IS REDUCED TO UPPER
C     TRAPEZOIDAL FORM F=QU BY IR ORTHOGONAL TRANSFORMATIONS
C     OF THE TYPE I-V*V**T/H WHERE **T DENOTES TRANSPOSE AND
C     F IS A WITH ITS COLUMNS PERMUTED. U IS GIVEN IN THE IR*N
C     UPPER TRAPEZIUM OF ARRAY A(IA,N). THE NON-ZERO ELEMENTS
C     OF THE VECTORS V, REQUIRED FOR THE FORMATION OF Q, ARE
C     GIVEN BELOW THE DIAGONAL IN THE FIRST IR COLUMNS OF
C     ARRAY A AND IN THE FIRST IR POSITIONS OF ARRAY D(N).
C     THE LOWER RIGHT-HAND SUBMATRIX OF A, OF ORDER
C     (M-IR)*(N-IR), IS OCCUPIED BY ELEMENTS WHICH CAN BE
C     REGARDED AS ZERO USING THE CRITERION TOL. AT EACH
C     STAGE OF THE REDUCTION THE ELEMENT OF LARGEST MODULUS
C     IN THE CURRENT MATRIX IS GIVEN IN AIJMAX AND INC RECORDS
C     THE COLUMN INTERCHANGE IN THE PIVOTING. A FAILURE EXIT
C     WILL BE DUE TO M LESS THAN N OR AN INVALID T.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01BKF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IA, IFAIL, IRANK, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), AIJMAX(N), D(N)
      INTEGER           INC(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  F, G, H, P, Q, SUMSQ, Z, ZERO
      INTEGER           I, ISAVE, J, K, KP1, L
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Data statements ..
      DATA              ZERO/0.0D0/
C     .. Executable Statements ..
C     NAG COPYRIGHT 1976
C     MARK 5 RELEASE
      ISAVE = IFAIL
      IFAIL = 1
      IF (M.GE.N) GO TO 20
      IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
   20 IFAIL = 2
      IF (T.GE.ZERO) GO TO 40
      IRANK = -1
      IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
C
C     QU DECOMPOSITION OF A WITH DETERMINATION OF RANK
C
   40 P = ZERO
      Z = ZERO
      DO 80 J = 1, N
         SUMSQ = ZERO
         DO 60 I = 1, M
            Q = ABS(A(I,J))
            SUMSQ = SUMSQ + Q**2
            IF (Q.GT.Z) Z = Q
   60    CONTINUE
         IF (P.GE.SUMSQ) GO TO 80
         P = SUMSQ
         L = J
   80 CONTINUE
      IF (Z.GT.T) GO TO 100
      IRANK = 0
      IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
  100 DO 260 K = 1, N
         AIJMAX(K) = Z
         INC(K) = L
         IF (P.NE.ZERO) GO TO 120
         IRANK = K - 1
         IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
         RETURN
C
C        INTERCHANGE COLUMNS K AND L OF A IF NECESSARY
C
  120    IF (L.EQ.K) GO TO 160
         DO 140 I = 1, M
            F = A(I,L)
            A(I,L) = A(I,K)
            A(I,K) = F
  140    CONTINUE
  160    F = A(K,K)
         H = P
         G = SQRT(H)
         IF (F.GE.ZERO) G = -G
         A(K,K) = G
         H = H - F*G
         F = F - G
         D(K) = F
         P = ZERO
         Z = ZERO
         IF (K.EQ.N) GO TO 240
         KP1 = K + 1
         DO 220 J = KP1, N
C
C           FORM JTH ELEMENT OF V**T*A
C
            G = F*A(K,J)
            DO 180 I = KP1, M
               G = G + A(I,K)*A(I,J)
  180       CONTINUE
            G = G/H
C
C           TRANSFORM JTH COLUMN OF A
C
            A(K,J) = A(K,J) - F*G
            SUMSQ = ZERO
            DO 200 I = KP1, M
               Q = A(I,J) - G*A(I,K)
               A(I,J) = Q
               Q = ABS(Q)
               IF (Q.GT.Z) Z = Q
               SUMSQ = SUMSQ + Q*Q
  200       CONTINUE
            IF (P.GE.SUMSQ) GO TO 220
            P = SUMSQ
            L = J
  220    CONTINUE
  240    IF (Z.LE.T) GO TO 280
  260 CONTINUE
C
C     RANK OF A IN IRANK
C
  280 IRANK = K
      IF (K.NE.N) AIJMAX(K+1) = Z
      IFAIL = 0
      RETURN
      END
