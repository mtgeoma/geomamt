      SUBROUTINE D05BYY(WT,IS,LENP,LIQ,SW,LSW,GAMFUN,AUXARY,A)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     --------------------------------------------------------------
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This routine evaluates the convolution term appearing in the
C     right hand side of the Vandermonde system by the FFT
C     techniques. these values are stored in SW.
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     --------------------------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER           IS, LENP, LIQ, LSW
C     .. Array Arguments ..
      DOUBLE PRECISION  A(2*LIQ), AUXARY(2*LIQ), GAMFUN(IS),
     *                  SW(LSW,0:IS-1), WT(LIQ)
C     .. Local Scalars ..
      DOUBLE PRECISION  RNUM, ST
      INTEGER           I, IFAIL, ISEC, M, NOWGT, NOWGT1, NOWGT2, S1
C     .. External Subroutines ..
      EXTERNAL          C06FAF, D05BYW, D05BYX
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
      IFAIL = 0
      NOWGT = LIQ
      NOWGT2 = 2*NOWGT
      NOWGT1 = NOWGT + 1
      S1 = IS - 1
C
      DO 20 I = 1, NOWGT
         AUXARY(I) = WT(I)
   20 CONTINUE
C
      DO 40 I = NOWGT1, NOWGT2
         AUXARY(I) = 0.0D0
   40 CONTINUE
C
      CALL C06FAF(AUXARY,NOWGT2,A,IFAIL)
C
      ST = 0.0D0
C
      DO 60 I = 1, NOWGT
         ST = ST + WT(I)
         SW(I+IS,0) = ST
   60 CONTINUE
C
      DO 140 M = 1, S1
         RNUM = 0.5D0*DBLE(M)
C
         DO 80 I = 1, NOWGT
            A(I) = DBLE((I+S1))**RNUM
   80    CONTINUE
C
         DO 100 I = NOWGT1, NOWGT2
            A(I) = 0.0D0
  100    CONTINUE
C
         ISEC = 0
         CALL D05BYW(A,AUXARY,NOWGT2,ISEC)
C
         DO 120 I = 1, NOWGT
            SW(I+IS,M) = A(I)
  120    CONTINUE
  140 CONTINUE
C
      CALL D05BYX(GAMFUN,IS)
C
      RETURN
      END
