      SUBROUTINE E04VDX(FIRSTV,NEGSTP,BIGALF,BIGBND,PNORM,JADD1,JADD2,
     *                  PALFA1,PALFA2,ISTATE,N,NCLIN0,NROWA,NCTOTL,
     *                  ANORM,AP,AX,BL,BU,FEATOL,P,X)
C     MARK 12 RE-ISSUE. NAG COPYRIGHT 1986.
C
C     ******************************************************************
C     E04VDX  FINDS STEPS  PALFA1, PALFA2  SUCH THAT
C     THE POINT  X + PALFA1*P  REACHES A LINEAR CONSTRAINT THAT IS
C                              CURRENTLY NOT IN THE WORKING SET BUT IS
C                              SATISFIED,
C     THE POINT  X + PALFA2*P  REACHES A LINEAR CONSTRAINT THAT IS
C                              CURRENTLY NOT IN THE WORKING SET BUT IS
C                              VIOLATED.
C     THE CONSTRAINTS ARE PERTURBED BY AN AMOUNT  FEATOL, SO THAT
C     PALFA1  IS SLIGHTLY LARGER THAN IT SHOULD BE, AND
C     PALFA2  IS SLIGHTLY SMALLER THAN IT SHOULD BE.  THIS GIVES
C     SOME LEEWAY LATER WHEN THE EXACT STEPS ARE COMPUTED BY E04VDW.
C
C     CONSTRAINTS IN THE WORKING SET ARE IGNORED  (ISTATE(J) GE 1).
C
C     IF  NEGSTP  IS TRUE, THE SEARCH DIRECTION WILL BE TAKEN TO BE
C     - P.
C
C
C     VALUES OF ISTATE(J)....
C
C     - 2         - 1         0           1          2         3
C     A*X LT BL   A*X GT BU   A*X FREE   A*X = BL   A*X = BU   BL = BU
C
C     THE VALUES -2 AND -1 DO NOT OCCUR ONCE E04MBY FINDS A FEASIBLE
C     POINT.
C
C     SYSTEMS OPTIMIZATION LABORATORY, STANFORD UNIVERSITY.
C     VERSION OF MAY 1982.  REV. OCT. 1982.
C     ******************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGALF, BIGBND, PALFA1, PALFA2, PNORM
      INTEGER           JADD1, JADD2, N, NCLIN0, NCTOTL, NROWA
      LOGICAL           FIRSTV, NEGSTP
C     .. Array Arguments ..
      DOUBLE PRECISION  ANORM(NCLIN0), AP(NCLIN0), AX(NROWA),
     *                  BL(NCTOTL), BU(NCTOTL), FEATOL(NCTOTL), P(N),
     *                  X(N)
      INTEGER           ISTATE(NCTOTL)
C     .. Scalars in Common ..
      INTEGER           ISTART, MSG, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  ABSATP, ATP, ATX, EPSPT9, ONE, RES, ROWNRM, ZERO
      INTEGER           I, J, JS
      LOGICAL           LASTV, NOLOW, NOUPP
C     .. Local Arrays ..
      CHARACTER*100     REC(3)
C     .. External Subroutines ..
      EXTERNAL          X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AE04VC/NOUT, MSG, ISTART
      COMMON            /CE04VC/PARM
C     .. Arrays in Common ..
      DOUBLE PRECISION  PARM(10)
C     .. Data statements ..
      DATA              ZERO/0.0D+0/, ONE/1.0D+0/
C     .. Executable Statements ..
C
      EPSPT9 = PARM(4)
      IF (MSG.EQ.99) THEN
         WRITE (REC,FMT=99999)
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         CALL X04BAF(NOUT,REC(3))
      END IF
      LASTV = .NOT. FIRSTV
      JADD1 = 0
      JADD2 = 0
      PALFA1 = BIGALF
      PALFA2 = ZERO
      IF (FIRSTV) PALFA2 = BIGALF
C
      DO 160 J = 1, NCTOTL
         JS = ISTATE(J)
         IF (JS.GT.0) GO TO 160
         IF (J.GT.N) GO TO 20
C
C        BOUND CONSTRAINT.
C
         ATX = X(J)
         ATP = P(J)
         ROWNRM = ONE
         GO TO 40
C
C        GENERAL LINEAR CONSTRAINT.
C
   20    I = J - N
         ATX = AX(I)
         ATP = AP(I)
         ROWNRM = ONE + ANORM(I)
C
   40    IF (NEGSTP) ATP = -ATP
         IF (ABS(ATP).GT.EPSPT9*ROWNRM*PNORM) GO TO 60
         RES = -ONE
         GO TO 140
   60    IF (ATP.GT.ZERO) GO TO 100
C
C        AX  IS DECREASING.
C        TEST FOR SMALLER PALFA1 IF LOWER BOUND IS SATISFIED.
C
         IF (JS.EQ.(-2)) GO TO 140
         ABSATP = -ATP
         NOLOW = BL(J) .LE. (-BIGBND)
         IF (NOLOW) GO TO 80
         RES = ATX - BL(J) + FEATOL(J)
         IF (BIGALF*ABSATP.LE.ABS(RES)) GO TO 80
         IF (PALFA1*ABSATP.LE.RES) GO TO 80
         PALFA1 = RES/ABSATP
         JADD1 = J
C
C        TEST FOR DIFFERENT PALFA2 IF UPPER BOUND IS VIOLATED.
C
   80    IF (JS.NE.(-1)) GO TO 140
         RES = ATX - BU(J) - FEATOL(J)
         IF (BIGALF*ABSATP.LE.ABS(RES)) GO TO 140
         IF (LASTV .AND. PALFA2*ABSATP.GE.RES) GO TO 140
         IF (FIRSTV .AND. PALFA2*ABSATP.LE.RES) GO TO 140
         PALFA2 = RES/ABSATP
         JADD2 = J
         GO TO 140
C
C        AX  IS INCREASING.
C        TEST FOR SMALLER PALFA1 IF UPPER BOUND IS SATISFIED.
C
  100    IF (JS.EQ.(-1)) GO TO 140
         NOUPP = BU(J) .GE. BIGBND
         IF (NOUPP) GO TO 120
         RES = BU(J) - ATX + FEATOL(J)
         IF (BIGALF*ATP.LE.ABS(RES)) GO TO 120
         IF (PALFA1*ATP.LE.RES) GO TO 120
         PALFA1 = RES/ATP
         JADD1 = J
C
C        TEST FOR DIFFERENT PALFA2 IF LOWER BOUND IS VIOLATED.
C
  120    IF (JS.NE.(-2)) GO TO 140
         RES = BL(J) - ATX - FEATOL(J)
         IF (BIGALF*ATP.LE.ABS(RES)) GO TO 140
         IF (LASTV .AND. PALFA2*ATP.GE.RES) GO TO 140
         IF (FIRSTV .AND. PALFA2*ATP.LE.RES) GO TO 140
         PALFA2 = RES/ATP
         JADD2 = J
C
  140    IF (MSG.EQ.99) THEN
            WRITE (REC,FMT=99998) J, JS, FEATOL(J), ATX, ATP, JADD1,
     *        PALFA1, JADD2, PALFA2
            CALL X04BAF(NOUT,REC(1))
         END IF
  160 CONTINUE
      RETURN
C
C
C     END OF E04VDX  ( BDPERT )
99999 FORMAT (/'    J  JS         FEATOL         AX             AP    ',
     *  ' JADD1       PALFA1     JADD2       PALFA2',/)
99998 FORMAT (I5,I4,3G15.5,2(I6,G17.7))
      END
