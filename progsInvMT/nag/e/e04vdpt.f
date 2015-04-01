      SUBROUTINE E04VDP(ORTHOG,UNITQ,INFORM,K1,K2,NACTIV,NCOLZ,NFREE,N,
     *                  NCTOTL,NQ,NROWA,NROWRT,NCOLRT,ISTATE,KACTIV,
     *                  KFREE,CONDMX,A,QTG,RT,ZY,WRK1,WRK2)
C     MARK 12 RE-ISSUE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-614 (APR 1988).
C
C *********************************************************************
C     E04VDP INCLUDES GENERAL LINEAR CONSTRAINTS  K1  THRU  K2  AS NEW
C     COLUMNS OF THE TQ FACTORIZATION STORED IN  RT, ZY.
C
C     SYSTEMS OPTIMIZATION LABORATORY, STANFORD UNIVERSITY.
C     VERSION OF SEPTEMBER 1981.  REV. OCT 1982, JAN 1983.
C *********************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CONDMX
      INTEGER           INFORM, K1, K2, N, NACTIV, NCOLRT, NCOLZ,
     *                  NCTOTL, NFREE, NQ, NROWA, NROWRT
      LOGICAL           ORTHOG, UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NROWA,N), QTG(N), RT(NROWRT,NCOLRT), WRK1(N),
     *                  WRK2(N), ZY(NQ,NQ)
      INTEGER           ISTATE(NCTOTL), KACTIV(N), KFREE(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, DTMAX, DTMIN
C     .. Local Scalars ..
      DOUBLE PRECISION  CSLAST, SNLAST
      INTEGER           I, IADD, IFX, ISWAP, JADD, K, L
C     .. External Subroutines ..
      EXTERNAL          E04VDZ
C     .. Common blocks ..
      COMMON            /HE04VC/ASIZE, DTMAX, DTMIN
C     .. Executable Statements ..
      DO 40 K = K1, K2
         IADD = KACTIV(K)
         JADD = N + IADD
         IF (NACTIV.EQ.NFREE) GO TO 20
C
         CALL E04VDZ(.FALSE.,.FALSE.,ORTHOG,UNITQ,INFORM,IFX,IADD,JADD,
     *               NACTIV,NCOLZ,NCOLZ,NFREE,N,NQ,NROWA,NROWRT,NCOLRT,
     *               KFREE,CONDMX,CSLAST,SNLAST,A,QTG,RT,ZY,WRK1,WRK2)
C
         IF (INFORM.GT.0) GO TO 20
         NACTIV = NACTIV + 1
         NCOLZ = NCOLZ - 1
         GO TO 40
C
   20    ISTATE(JADD) = 0
         KACTIV(K) = -KACTIV(K)
   40 CONTINUE
C
      IF (NACTIV.EQ.K2) RETURN
C
C     SOME OF THE CONSTRAINTS WERE CLASSED AS DEPENDENT AND NOT INCLUDED
C     IN THE FACTORIZATION.  MOVE ACCEPTED INDICES TO THE FRONT OF
C     KACTIV  AND SHIFT REJECTED INDICES (WITH NEGATIVE VALUES) TO
C     THE END.
C
      L = K1 - 1
      DO 60 K = K1, K2
         I = KACTIV(K)
         IF (I.LT.0) GO TO 60
         L = L + 1
         IF (L.EQ.K) GO TO 60
         ISWAP = KACTIV(L)
         KACTIV(L) = I
         KACTIV(K) = ISWAP
   60 CONTINUE
      RETURN
C
C     END OF E04VDP (TQADD)
      END
