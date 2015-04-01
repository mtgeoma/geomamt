      SUBROUTINE F01BQZ(N,EPS,RL,IRL,D,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     This routine originally called F01BQF.
C     FORMCHOL AND FIND BOUND
C
C     GIVEN AN N*N SYMMETRIC MATRIX G THIS SUBROUTINE FORMS THE
C     CHOLESKY FACTORIZATION L*D*LT OF A POSITIVE-DEFINITE
C     MATRIX G + E, WHERE LT IS THE TRANSPOSE OF L. L IS AN N*N
C     UNIT LOWER-TRIANGULAR MATRIX STORED ROW BY ROW OMITTING
C     THE UNIT DIAGONAL, AND D IS A DIAGONAL MATRIX. THE
C     MATRIX E IS DIAGONAL AND IS IDENTICALLY ZERO WHEN G IS
C     SUFFICIENTLY POSITIVE DEFINITE. THE MATRIX G IS STORED
C     WITH ITS DIAGONAL IN THE ARRAY D(I), I = 1(1)N, AND
C     ITS LOWER TRIANGLE STORED ROW BY ROW IN THE ARRAY RL(I)
C     I = 1(1)N(N-1)/2. THE FACTORIZATION IS OVERWRITTEN ON
C     RL AND D. EPS SHOULD BE SET TO EPSMCH AT ENTRY
C     AND ON EXIT WILL CONTAIN EPSMCH*THE LARGER OF 1.0 AND THE
C     ELEMENT OF G WHICH IS LARGEST IN MODULUS. IF E IS NOT A ZERO
C     MATRIX THEN IFAIL CORRESPONDS TO THE POSITION OF THE MOST
C     NEGATIVE
C     ELEMENT OF THE DIAGONAL MATRIX D - E. IF THE MATRIX E IS ZERO
C     THEN IFAIL = 0.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, SHEDMOND R.
C     GRAHAM
C     AND MARGARET H. WRIGHT
C     D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND.
C     CROWN COPYRIGHT RESERVED
C     SEPTEMBER 1975
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01BQZ')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS
      INTEGER           IFAIL, IRL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(N), RL(IRL)
C     .. Local Scalars ..
      DOUBLE PRECISION  BOUND, DEL, EPSMCH, GAMMA, ONE, RMAX, U, V, W,
     *                  WMU, XI, XIN, Z, ZERO
      INTEGER           I, IQ, IR, IS, ISAVE, IT, J, JP1, K, M, NH
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Data statements ..
      DATA              ZERO/0.0D+0/, ONE/1.0D+0/
C     .. Executable Statements ..
      ISAVE = IFAIL
      EPSMCH = EPS
C
C     FIND BOUND, WHICH IS THE MAXIMUM OF ABS(RL(I))/N
C     ABS(D(I)) AND EPSMCH.
C
      IF (N.LE.1) GO TO 40
      XI = ABS(RL(1))
      IF (N.LE.2) GO TO 60
      NH = N*(N-1)/2
      DO 20 I = 2, NH
         U = ABS(RL(I))
         IF (U.GT.XI) XI = U
   20 CONTINUE
      GO TO 60
   40 XI = ZERO
   60 XIN = N
      XIN = XI/XIN
      GAMMA = ABS(D(1))
      IF (N.EQ.1) GO TO 100
      DO 80 I = 2, N
         U = ABS(D(I))
         IF (U.GT.GAMMA) GAMMA = U
   80 CONTINUE
  100 IF (XIN.LE.GAMMA) GO TO 120
      BOUND = XIN
      GO TO 140
  120 BOUND = GAMMA
  140 CONTINUE
      IF (BOUND.LT.EPSMCH) BOUND = EPSMCH
      IF (XI.LE.GAMMA) GO TO 160
      DEL = XI*EPSMCH
      GO TO 180
  160 DEL = GAMMA*EPSMCH
  180 IF (DEL.LT.EPSMCH) DEL = EPSMCH
      IFAIL = 0
      M = 1
      Z = ONE/BOUND
      RMAX = ZERO
      DO 360 J = 1, N
         U = D(J)
         IQ = M
         IF (J.EQ.1) GO TO 220
         IT = J - 1
         DO 200 K = 1, IT
            V = RL(M)
            W = V/D(K)
            RL(M) = W
            M = M + 1
            U = U - W*V
  200    CONTINUE
  220    IF (U.GE.DEL) GO TO 240
         W = ABS(U)
         IF (W.LE.DEL) W = DEL
         WMU = W - U
         IF (WMU.LT.RMAX) GO TO 260
         RMAX = WMU
         IFAIL = J
         GO TO 260
  240    W = U
  260    V = ZERO
         IS = M
         IF (J.EQ.N) GO TO 340
         JP1 = J + 1
         DO 320 I = JP1, N
            U = ZERO
            IR = IQ
            IF (J.EQ.1) GO TO 300
            IT = J - 1
            DO 280 K = 1, IT
               U = U - RL(IS)*RL(IR)
               IS = IS + 1
               IR = IR + 1
  280       CONTINUE
  300       RL(IS) = U + RL(IS)
            U = RL(IS)
            IS = IS + I - J
            U = ABS(U)
            IF (U.GT.V) V = U
  320    CONTINUE
C
C        D(J) ARE MODIFIED SUCH THAT ABS( SQRT(D(J)) * ANY ELEMENT
C        OF THE JTH COLUMN OF L ) LE SQRT( BOUND )
C
  340    V = V*V*Z
         D(J) = W
         IF (W.LT.V) D(J) = V
  360 CONTINUE
      EPS = DEL
      IF (IFAIL.NE.0) IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
C     END OF SUBROUTINE F01BQZ
C
      END
