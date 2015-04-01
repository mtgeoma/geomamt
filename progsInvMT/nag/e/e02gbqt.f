      SUBROUTINE E02GBQ(K,N,ZZ,IZR,IRR,DD,RR,IC)
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
C     AND GIVEN THE INDEX  IC  OF A COLUMN
C     TO BE REMOVED  (1 .LE. IC .LE. K),
C     THIS ROUTINE UPDATES THE FACTORIZATION.
C
C     THE VALUE OF  K  IS DECREASED BY ONE.
C     ***************
C
C     .. Scalar Arguments ..
      INTEGER           IC, IRR, IZR, K, N
C     .. Array Arguments ..
      DOUBLE PRECISION  DD(N), RR(IRR), ZZ(IZR,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  DI, RJ, SRA, SRB, ST1, ST2
      INTEGER           I, IC1, IF1, IF2, IM1, J, JEND, JINC, JSTRT,
     *                  KM1, LSTRT
C     .. External Subroutines ..
      EXTERNAL          E02GBN, E02GBP
C     .. Executable Statements ..
      KM1 = K - 1
      IC1 = IC + 1
C
C     ***************
C     SPECIAL CASES ARE HANDLED FIRST.
C     1.  K=1 AND THE FACTORIZATION BECOMES NULL.
C     2.  IC=K AND THE UPDATING IS TRIVIAL.
C     ***************
C
      IF (K.GT.1) GO TO 20
      K = 0
      RETURN
   20 CONTINUE
      IF (IC.LT.K) GO TO 40
      K = KM1
      RETURN
   40 CONTINUE
C
C     ***************
C     GENERAL UPDATING STEP.
C     THE COLUMN TO BE DELETED MUST BE PERMUTED
C     TO THE RIGHT, AND SUBDIAGONAL ELEMENTS
C     WHICH RESULT IN  RR  HAVE TO BE
C     TRANSFORMED TO ZERO.
C     ***************
C
      JSTRT = IC1
      JEND = K
      JINC = N
      DO 100 I = 1, K
         LSTRT = JSTRT - JINC - 1
         DO 60 J = JSTRT, JEND
            RR(J) = RR(J+1)
   60    CONTINUE
         IF (I.LE.IC) GO TO 80
         IM1 = I - 1
         DI = DD(I)
         RJ = RR(JSTRT)
         CALL E02GBP(DD(IM1),DI,RR(LSTRT),RJ,ST1,ST2,SRA,SRB,IF1,IF2)
         RR(JSTRT) = RJ
         DD(I) = DI
         IF (JEND.GT.JSTRT) CALL E02GBN(JEND-JSTRT,RR(LSTRT+1)
     *                                  ,1,RR(JSTRT+1),1,ST1,ST2,SRA,
     *                                  SRB,IF1,IF2)
         CALL E02GBN(N,ZZ(1,IM1),1,ZZ(1,I),1,ST1,ST2,SRA,SRB,IF1,IF2)
   80    CONTINUE
         JSTRT = JSTRT + JINC
         JEND = JEND + JINC
         JINC = JINC - 1
         IF (I.GT.IC) JSTRT = JSTRT + 1
  100 CONTINUE
      K = KM1
      RETURN
      END
