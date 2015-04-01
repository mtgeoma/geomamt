      SUBROUTINE D02QZZ(T,YINT,YPINT,NINT,NEQ,X,YY,P,PHI)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C
C     MODIFIED VERSION OF D02QFS - FIRST NINT COMPONENTS INTERPOLATED
C                                  ----------
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, X
      INTEGER           NEQ, NINT
C     .. Array Arguments ..
      DOUBLE PRECISION  P(NEQ), PHI(NEQ,*), YINT(NINT), YPINT(NINT),
     *                  YY(NEQ)
C     .. Scalars in Common ..
      DOUBLE PRECISION  HOLD, TOLD, TSTAR, TWOU, XOLD, XSAVE
      INTEGER           IBEGIN, IINTEG, INFLOP, INIT, IQUIT, ITOL,
     *                  ITSTOP, IVC, KGI, KLE4, KOLD, KORD, KPREV,
     *                  KSTEPS, NS
      LOGICAL           INTOUT, NORND, PHASE1, START, STIFF
C     .. Arrays in Common ..
      DOUBLE PRECISION  ALPHA(12), BETA(12), G(13), GI(11), PSI(12),
     *                  SIG(13), V(12), W(12)
      INTEGER           IV(10)
C     .. External Subroutines ..
      EXTERNAL          D02QFR
C     .. Common blocks ..
      COMMON            /AD02QF/ALPHA, BETA, PSI, V, W, SIG, G, GI,
     *                  XOLD, HOLD, TOLD, XSAVE, TSTAR, TWOU, INIT,
     *                  IBEGIN, ITOL, IINTEG, ITSTOP, INFLOP, IQUIT, IV,
     *                  NS, KORD, KOLD, KSTEPS, KLE4, KPREV, IVC, KGI,
     *                  START, PHASE1, NORND, STIFF, INTOUT
C     .. Save statement ..
      SAVE              /AD02QF/
C     .. Executable Statements ..
C
      CALL D02QFR(X,YY,T,YINT,YPINT,NINT,NEQ,KOLD,PHI,IVC,IV,KGI,GI,
     *            ALPHA,G,W,XOLD,P)
C
      RETURN
      END
