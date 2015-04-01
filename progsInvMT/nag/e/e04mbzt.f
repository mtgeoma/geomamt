      SUBROUTINE E04MBZ(N,NACTIV,NCTOTL,NFREE,JBIGST,KBIGST,ISTATE,
     *                  KACTIV,DINKY,FEAMIN,TRULAM,FEATOL,RLAMDA)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C *********************************************************************
C     FIND THE BIGGEST SCALED MULTIPLIER LARGER THAN UNITY.
C
C     SYSTEMS OPTIMIZATION LABORATORY, STANFORD UNIVERSITY.
C     ORIGINAL VERSION DECEMBER 1982.
C *********************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DINKY, FEAMIN, TRULAM
      INTEGER           JBIGST, KBIGST, N, NACTIV, NCTOTL, NFREE
C     .. Array Arguments ..
      DOUBLE PRECISION  FEATOL(NCTOTL), RLAMDA(N)
      INTEGER           ISTATE(NCTOTL), KACTIV(N)
C     .. Scalars in Common ..
      INTEGER           ISTART, MSG, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGGST, ONE, RLAM
      INTEGER           IS, J, K, NFIXED, NLAM
C     .. Local Arrays ..
      CHARACTER*33      REC(3)
C     .. External Subroutines ..
      EXTERNAL          X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AE04VC/NOUT, MSG, ISTART
C     .. Data statements ..
      DATA              ONE/1.0D+0/
C     .. Executable Statements ..
C
      JBIGST = 0
      NFIXED = N - NFREE
      NLAM = NFIXED + NACTIV
      IF (NLAM.EQ.0) GO TO 40
C
      BIGGST = ONE + DINKY
      DO 20 K = 1, NLAM
         J = KACTIV(K)
         IF (K.LE.NACTIV) J = J + N
         IS = ISTATE(J)
         IF (IS.LT.1) GO TO 20
         RLAM = RLAMDA(K)
         IF (IS.EQ.2) RLAM = -RLAM
         IF (IS.EQ.3) RLAM = ABS(RLAM)
         RLAM = (FEATOL(J)/FEAMIN)*RLAM
C
         IF (BIGGST.GE.RLAM) GO TO 20
         BIGGST = RLAM
         TRULAM = RLAMDA(K)
         JBIGST = J
         KBIGST = K
   20 CONTINUE
      IF (MSG.GE.80) THEN
         WRITE (REC,FMT=99999) JBIGST, BIGGST
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         CALL X04BAF(NOUT,REC(3))
      END IF
C
   40 RETURN
C
C
C     END OF E04MBZ  ( LPBGST )
99999 FORMAT (/' //E04MBZ// JBIGST         BIGGST',/' //E04MBZ//  ',I5,
     *  G15.4)
      END
