      SUBROUTINE E04MFK(MSGLVL,N,NCLIN,NCTOTL,BIGBND,NAMED,NAMES,ISTATE,
     *                  BL,BU,CLAMDA,FEATOL,R)
C     MARK 17 RE-ISSUE. NAG COPYRIGHT 1995.
C
C     ==================================================================
C     E04MFK  prints r(x) (x, A*x and c(x)), the bounds, the
C     multipliers, and the slacks (distance to the nearer bound).
C
C     Original Fortran 77 version written  October 1984.
C     This version of  E04MFK dated  05-May-93.
C     ==================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND
      INTEGER           MSGLVL, N, NCLIN, NCTOTL
      LOGICAL           NAMED
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(NCTOTL), BU(NCTOTL), CLAMDA(NCTOTL),
     *                  FEATOL(NCTOTL), R(NCTOTL)
      INTEGER           ISTATE(NCTOTL)
      CHARACTER*8       NAMES(*)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  B1, B2, RJ, SLK, SLK1, SLK2, TOL, WLAM
      INTEGER           IS, J, NPLIN, NUMBER
      CHARACTER         KEY
      CHARACTER*2       STATE
      CHARACTER*8       NAME
      CHARACTER*80      LINE
C     .. Local Arrays ..
      CHARACTER*2       LSTATE(-2:4)
      CHARACTER*80      REC(4)
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Data statements ..
      DATA              LSTATE(-2)/'--'/, LSTATE(-1)/'++'/
      DATA              LSTATE(0)/'FR'/, LSTATE(1)/'LL'/
      DATA              LSTATE(2)/'UL'/, LSTATE(3)/'EQ'/
      DATA              LSTATE(4)/'TF'/
C     .. Executable Statements ..
C
      IF (MSGLVL.LT.10 .AND. MSGLVL.NE.1) RETURN
C
      WRITE (REC,FMT=99999) 'Varbl'
      CALL X04BAY(IPRINT,4,REC)
      NAME = 'V'
      NPLIN = N + NCLIN
C
      DO 20 J = 1, NCTOTL
         B1 = BL(J)
         B2 = BU(J)
         WLAM = CLAMDA(J)
         RJ = R(J)
C
         IF (J.LE.N) THEN
            NUMBER = J
         ELSE IF (J.LE.NPLIN) THEN
            NUMBER = J - N
            IF (NUMBER.EQ.1) THEN
               WRITE (REC,FMT=99999) 'L Con'
               CALL X04BAY(IPRINT,4,REC)
               NAME = 'L'
            END IF
         ELSE
            NUMBER = J - NPLIN
            IF (NUMBER.EQ.1) THEN
               WRITE (REC,FMT=99999) 'N Con'
               CALL X04BAY(IPRINT,4,REC)
               NAME = 'N'
            END IF
         END IF
C
C        Print a line for the jth variable or constraint.
C        ------------------------------------------------
         IS = ISTATE(J)
         STATE = LSTATE(IS)
         TOL = FEATOL(J)
         SLK1 = RJ - B1
         SLK2 = B2 - RJ
         IF (ABS(SLK1).LT.ABS(SLK2)) THEN
            SLK = SLK1
            IF (B1.LE.-BIGBND) SLK = SLK2
         ELSE
            SLK = SLK2
            IF (B2.GE.BIGBND) SLK = SLK1
         END IF
C
C        Flag infeasibilities, primal and dual degeneracies,
C        and active QP constraints that are loose in NP.
C
         KEY = ' '
         IF (SLK1.LT.-TOL .OR. SLK2.LT.-TOL) KEY = 'I'
         IF (IS.EQ.0 .AND. ABS(SLK).LE.TOL) KEY = 'D'
         IF (IS.GE.1 .AND. ABS(WLAM).LE.TOL) KEY = 'A'
C
         WRITE (LINE,FMT=99998) NAME, NUMBER, KEY, STATE, RJ, B1, B2,
     *     WLAM, SLK
C
C        Reset special cases:
C           Infinite bounds
C           Zero bounds
C           Lagrange multipliers for inactive constraints
C           Lagrange multipliers for infinite bounds
C           Infinite slacks
C           Zero slacks
C
         IF (B1.LE.-BIGBND) LINE(28:41) = '      None   '
         IF (B2.GE.BIGBND) LINE(42:55) = '      None   '
         IF (B1.EQ.ZERO) LINE(28:41) = '       .     '
         IF (B2.EQ.ZERO) LINE(42:55) = '       .     '
         IF (IS.EQ.0 .OR. WLAM.EQ.ZERO) LINE(56:67) = '       .   '
         IF (B1.LE.-BIGBND .AND. B2.GE.BIGBND) THEN
            LINE(56:67) = '           '
            LINE(68:79) = '           '
         END IF
         IF (SLK.EQ.ZERO) LINE(68:79) = '     .     '
         CALL X04BAF(IPRINT,LINE)
C
   20 CONTINUE
C
      RETURN
C
C
C     End of E04MFK.  (CMPRNT)
C
99999 FORMAT (//1X,A5,1X,'State',5X,'Value',7X,'Lower Bound',3X,'Upper',
     *       ' Bound',4X,'Lagr Mult',3X,'Slack',/)
99998 FORMAT (1X,A1,I4,2X,A1,1X,A2,1X,1P,3G14.6,1P,2G12.4)
      END
