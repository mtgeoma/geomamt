      SUBROUTINE F01LBF(N,M1,M2,A,IA,AL,IL,IN,IV,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 14 REVISED. IER-727 (DEC 1989).
C
C     PDK, DNAC, NPL, TEDDINGTON. JUNE 1976.
C     NPL DNAC LIBRARY SUBROUTINE BANDLU.
C     REVISED BY N.A.G. CENTRAL OFFICE, 1979.
C
C     LU DECOMPOSITION OF A BAND MATRIX A WITH M1 SUB AND M2 SUPER-
C     DIAGONALS. THE MATRIX A IS PRESENTED AS AN IW*N ARRAY
C     WITH EACH ROW TOP JUSTIFIED, WHERE IW=MIN(N,M1+M2+1).
C     L AND U ARE FOUND USING GAUSSIAN
C     ELIMINATION WITH PARTIAL PIVOTING. U IS OVERWRITTEN ON A, THE
C     DIAGONAL ELEMENTS BEING STORED AS THEIR RECIPROCALS AND L IS
C     STORED AS A SEPARATE N*M1 MATRIX AL. THE ARRAY AL MUST
C     HAVE DIMENSIONS AT LEAST N*1 IN THE CALLING PROGRAM
C     IF M1 IS ZERO. DETAILS OF THE PIVOTING ARE
C     STORED IN THE VECTOR IN WHICH IS SUCH THAT IN(I)=K IF ROWS I
C     AND K WERE INTERCHENGED AT THE ITH MAJOR STEP.
C     IA AND IL ARE THE ROW DIMENSIONS OF A AND AL IN THE
C     CALLING PROGRAM.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01LBF')
C     .. Scalar Arguments ..
      INTEGER           IA, IFAIL, IL, IV, M1, M2, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), AL(IL,*)
      INTEGER           IN(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, ONE, X, Y, ZERO
      INTEGER           I, IK, IR, ISAVE, IW, IWW, J, JR, K, M
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Data statements ..
      DATA              ZERO/0.0D0/, ONE/1.0D0/
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL = 1
      IV = 0
      EPS = X02AJF()
      IW = MIN(M1+M2+1,N)
      IF (N.LT.1 .OR. M1.LT.0 .OR. M1.GT.N-1 .OR. M2.LT.0 .OR. M2.GT.N-
     *    1 .OR. IA.LT.IW .OR. IL.LT.1 .OR. IL.LT.M1) GO TO 380
      IFAIL = 2
      M = M2 + 1
      K = IW - M2 - 1
      IF (K.LE.0) GO TO 60
      DO 40 I = 1, K
         M = M + 1
         DO 20 J = M, IW
            A(J,I) = ZERO
   20    CONTINUE
   40 CONTINUE
   60 M = N - IW + M1 + 2
      J = IW + 1
      IF (M.GT.N) GO TO 120
      DO 100 I = M, N
         J = J - 1
         DO 80 K = J, IW
            A(K,I) = ZERO
   80    CONTINUE
  100 CONTINUE
C
C     ZEROS INSERTED.
C
  120 DO 180 I = 1, N
         X = ZERO
         DO 140 J = 1, IW
            X = X + ABS(A(J,I))
  140    CONTINUE
         IF (X.GT.ZERO) GO TO 160
         IR = I
         GO TO 360
  160    AL(1,I) = ONE/X
  180 CONTINUE
C
C     ROW NORMS OF A CALCULATED AND THEIR RECIPROCALS
C     STORED IN FIRST COLUMN OF AL.
C
      IFAIL = 3
      DO 340 IR = 1, N
         X = ZERO
         M = MIN(N,IR+M1)
         IWW = MIN(IW,N-IR+1)
         DO 200 I = IR, M
            Y = ABS(A(1,I))*AL(1,I)
            IF (Y.LE.X) GO TO 200
            X = Y
            J = I
  200    CONTINUE
         IF (X.LT.EPS) GO TO 360
         IN(IR) = J
C
C        (IR)TH PIVOT ELEMENT SELECTED.
C
         IF (J.EQ.IR) GO TO 240
         DO 220 I = 1, IWW
            X = A(I,IR)
            A(I,IR) = A(I,J)
            A(I,J) = X
  220    CONTINUE
         AL(1,J) = AL(1,IR)
C
C        ROWS IR AND J INTERCHANGED.
C
  240    JR = IR + 1
         Y = ONE/A(1,IR)
         IF (JR.GT.M) GO TO 320
         DO 300 I = JR, M
            X = A(1,I)*Y
            IF (IWW.LT.2) GO TO 280
            DO 260 J = 2, IWW
               A(J-1,I) = A(J,I) - X*A(J,IR)
  260       CONTINUE
  280       IK = I - IR
            AL(IK,IR) = X
            A(IWW,I) = ZERO
  300    CONTINUE
  320    A(1,IR) = Y
  340 CONTINUE
C
C     ELIMINATION COMPLETED.
C
      IFAIL = 0
      RETURN
  360 A(1,IR) = ZERO
      IV = IR
  380 IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
      END
