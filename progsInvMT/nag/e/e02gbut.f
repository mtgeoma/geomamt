      SUBROUTINE E02GBU(K,NCOL,N,NUM,E,IER,MPL1,RES,P,PTA,INDX,ALF)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 8A REVISED. IER-256 (AUG 1980).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     ***************
C     THIS ROUTINE DETERMINES ALL OF THE RATIOS  ALF
C     OF THE FORM
C     -RES(I)/(E(.,I)-TRANSP)*P,
C     FOR  I = K+1,...,NCOL
C     WHICH ARE POSITIVE AND HENCE INDICATE DISTANCES
C     FROM THE POINT  X  TO BREAKPOINTS WHICH WILL
C     BE ENCOUNTERED IN TRAVEL ALONG DIRECTION  P.
C     THE INDEX VECTOR  INDX  IS REARRANGED SO THAT
C     ITS  K+1  THROUGH  NUM  COMPONENTS CORRESPOND TO
C     THESE POSITIVE RATIOS.
C     THE RESULTS ARE HEAPED SO THAT THE  ALF  VALUES CAN
C     BE INSPECTED IN ORDER FROM SMALLEST TO LARGEST.
C
C     THE INNER PRODUCTS (E(.,I)-TRANSPOSE)*P ARE SAVED
C     FOR LATER USE IN UPDATING THE RESIDUAL VALUES.
C
C     ***************
C
C     .. Scalar Arguments ..
      INTEGER           IER, K, MPL1, N, NCOL, NUM
C     .. Array Arguments ..
      DOUBLE PRECISION  ALF(NCOL), E(IER,MPL1), P(N), PTA(NCOL),
     *                  RES(NCOL)
      INTEGER           INDX(NCOL)
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPS
C     .. Local Scalars ..
      DOUBLE PRECISION  BIG, DEN, ONE, RATIO, RESID, ZERO
      INTEGER           I, IX, JX, KP1, NUMK
C     .. External Functions ..
      DOUBLE PRECISION  E02GBJ, X02ALF
      EXTERNAL          E02GBJ, X02ALF
C     .. External Subroutines ..
      EXTERNAL          E02GBM
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AE02GB/EPS
C     .. Data statements ..
      DATA              ZERO/0.0D+00/
      DATA              ONE/1.0D+00/
C     .. Executable Statements ..
      BIG = X02ALF()
      NUM = 0
      KP1 = K + 1
      IF (KP1.GT.NCOL) GO TO 60
      DO 40 I = KP1, NCOL
         IX = INDX(I)
         RESID = RES(IX)
         DEN = E02GBJ(N,E(1,IX),1,P,1,N,N)
         PTA(IX) = DEN
         IF (RESID*DEN.GE.ZERO) GO TO 40
         RESID = ABS(RESID)
         DEN = ABS(DEN)
         IF (DEN.GE.ONE) GO TO 20
         IF (RESID.GE.DEN*BIG) GO TO 40
   20    CONTINUE
         RATIO = RESID/DEN
         NUM = NUM + 1
         NUMK = NUM + K
         JX = INDX(NUMK)
         INDX(NUMK) = IX
         INDX(I) = JX
         ALF(NUM) = RATIO
   40 CONTINUE
C
C     ***************
C     HEAP THE POSITIVE RATIOS
C     ***************
C
      IF (NUM.GT.0) CALL E02GBM(.TRUE.,NUM,INDX(KP1),ALF,NUM)
   60 CONTINUE
      IF (1.GT.K) GO TO 100
      DO 80 I = 1, K
         IX = INDX(I)
         PTA(IX) = E02GBJ(N,E(1,IX),1,P,1,N,N)
   80 CONTINUE
  100 CONTINUE
      RETURN
      END
