      SUBROUTINE E02GBR(K,N,ZZ,IZR,IRR,DD,RR,COL,W,IW)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     ***************
C     GIVEN THE FACTORIZATION
C     ZZ*DD*RR
C     OF SOME  N BY K  MATRIX, WHERE
C     (ZZ-TRANSP)*(ZZ) = (DD-INV),
C     DD  IS DIAGONAL AND NONSINGULAR,
C     AND
C     RR  HAS ZEROS BELOW THE DIAGONAL,
C     AND GIVEN A  K+1ST  COLUMN TO BE ADDED,
C     THIS ROUTINE UPDATES THE FACTORIZATION.
C
C     W  IS A SCRATCH VECTOR.
C     ***************
C
C     .. Scalar Arguments ..
      INTEGER           IRR, IW, IZR, K, N
C     .. Array Arguments ..
      DOUBLE PRECISION  COL(N), DD(N), RR(IRR), W(IW), ZZ(IZR,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  DI, ONE, SRA, SRB, ST1, ST2, WI, ZERO
      INTEGER           I, I1, IF1, IF2, J, JDEL, KP1, KP2, NM1
C     .. External Functions ..
      DOUBLE PRECISION  E02GBJ
      EXTERNAL          E02GBJ
C     .. External Subroutines ..
      EXTERNAL          E02GBN, E02GBP
C     .. Data statements ..
      DATA              ZERO/0.0D+00/, ONE/1.0D+00/
C     .. Executable Statements ..
      IF (K.GT.0) GO TO 80
C
C     ***************
C     FOR THE SPECIAL CASE THAT   ZZ*DD*RR
C     WAS VACUOUS,  RESET THE ARRAYS  ZZ
C     AND  DD  TO IDENTITIES IN ORDER TO
C     BEGIN THE PROCESS.
C     ***************
C
      K = 0
      NM1 = N - 1
      IF (1.GT.NM1) GO TO 60
      DO 40 I = 1, NM1
         I1 = I + 1
         DO 20 J = I1, N
            ZZ(I,J) = ZERO
            ZZ(J,I) = ZERO
   20    CONTINUE
         ZZ(I,I) = ONE
         DD(I) = ONE
   40 CONTINUE
   60 CONTINUE
      ZZ(N,N) = ONE
      DD(N) = ONE
   80 CONTINUE
      KP1 = K + 1
      KP2 = K + 2
C
C     ***************
C     TRANSFORM THE INCOMING COLUMN,
C     AND STORE THE RESULT IN  W.
C     ***************
C
      DO 100 I = 1, N
         W(I) = E02GBJ(N,ZZ(1,I),1,COL,1,N,N)
  100 CONTINUE
C
C     ***************
C     ZERO OUT THE SPIKE WHICH RESULTS FROM
C     STORING  W  IN  RR.   UPDATE  ZZ  AND  DD.
C     ***************
C
      IF (KP2.GT.N) GO TO 140
      DO 120 I = KP2, N
         DI = DD(I)
         WI = W(I)
         CALL E02GBP(DD(KP1),DI,W(KP1),WI,ST1,ST2,SRA,SRB,IF1,IF2)
         W(I) = WI
         DD(I) = DI
         CALL E02GBN(N,ZZ(1,KP1),1,ZZ(1,I),1,ST1,ST2,SRA,SRB,IF1,IF2)
  120 CONTINUE
  140 CONTINUE
C
C     ***************
C     STORE THE RESULTING NEW COLUMN IN  RR.
C     ***************
C
      J = KP2
      JDEL = N
      DO 160 I = 1, KP1
         RR(J) = W(I)
         J = J + JDEL
         JDEL = JDEL - 1
  160 CONTINUE
      K = KP1
      RETURN
      END
