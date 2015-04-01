      SUBROUTINE H02BBW(N,NCTOTL,A,LDA,BL1,BU1,X,OBJMIP,IOPTCL,ISTATE,
     *                  CLAM,NODKNT,ISTAT,BIGBND)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-928 (APR 1991).
C     MARK 16 RE-ISSUE. NAG COPYRIGHT 1992.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, OBJMIP
      INTEGER           IOPTCL, ISTAT, LDA, N, NCTOTL, NODKNT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), BL1(NCTOTL), BU1(NCTOTL),
     *                  CLAM(NCTOTL), X(N)
      INTEGER           ISTATE(NCTOTL,IOPTCL)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  B1, B2, RES, RES2, V, WLAM
      INTEGER           IS, JK, KK
      CHARACTER*2       LS
C     .. Local Arrays ..
      CHARACTER*2       LSTATE(7)
      CHARACTER*80      REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
C     .. Data statements ..
      DATA              LSTATE(1)/'--'/, LSTATE(2)/'++'/
      DATA              LSTATE(3)/'FR'/, LSTATE(4)/'LL'/
      DATA              LSTATE(5)/'UL'/, LSTATE(6)/'EQ'/
      DATA              LSTATE(7)/'TB'/
C     .. Executable Statements ..
C
      IF (ISTAT.EQ.5) GO TO 20
      WRITE (REC,FMT=99985) ABS(NODKNT)
      CALL X04BAY(NOUT,4,REC)
      IF (ISTAT.LT.1) THEN
         WRITE (REC,FMT=99999)
      ELSE
         IF (ISTAT.EQ.1) WRITE (REC,FMT=99998)
         IF (ISTAT.GT.1) WRITE (REC,FMT=99997)
      END IF
C
      CALL X04BAY(NOUT,2,REC)
      WRITE (REC,FMT=99996) OBJMIP
      CALL X04BAY(NOUT,2,REC)
      IF (NODKNT.LT.0) THEN
         NODKNT = ABS(NODKNT)
         RETURN
      END IF
C
   20 WRITE (REC,FMT=99995)
      CALL X04BAY(NOUT,4,REC)
C
      DO 40 JK = 1, NCTOTL
         B1 = BL1(JK)
         B2 = BU1(JK)
         WLAM = CLAM(JK)
         IS = ISTATE(JK,IOPTCL)
         LS = LSTATE(IS+3)
         IF (JK.LE.N) THEN
C
C           Section 1 -- the variables  x.
C           ------------------------------
            KK = JK
            V = X(JK)
C
         ELSE
            IF (JK.EQ.N+1) THEN
               WRITE (REC,FMT=99994)
               CALL X04BAY(NOUT,4,REC)
            END IF
C
C
C           Section 2 -- the linear constraints  A*x.
C           -----------------------------------------
C
            KK = JK - N
            V = DDOT(N,A(KK,1),LDA,X,1)
C
         END IF
C
C        Print a line for the JK-th variable or constraint.
C        -------------------------------------------------
         RES = V - B1
         RES2 = B2 - V
         IF (ABS(RES).GT.ABS(RES2)) RES = RES2
C
         IF (JK.LE.N) THEN
            WRITE (REC,FMT=99993) KK, LS, V, B1, B2, WLAM, RES
            IF (B1.LE.-BIGBND) WRITE (REC,FMT=99992) KK, LS, V, B2,
     *          WLAM, RES
            IF (B2.GE.BIGBND) WRITE (REC,FMT=99991) KK, LS, V, B1, WLAM,
     *          RES
            IF (B1.LE.-BIGBND .AND. B2.GE.BIGBND) WRITE (REC,FMT=99990)
     *          KK, LS, V, WLAM, RES
            CALL X04BAF(NOUT,REC(1))
C
         ELSE
            WRITE (REC,FMT=99986) KK, LS, V, B1, B2, WLAM, RES
            IF (B1.LE.-BIGBND) WRITE (REC,FMT=99989) KK, LS, V, B2,
     *          WLAM, RES
            IF (B2.GE.BIGBND) WRITE (REC,FMT=99988) KK, LS, V, B1, WLAM,
     *          RES
            IF (B1.LE.-BIGBND .AND. B2.GE.BIGBND) WRITE (REC,FMT=99987)
     *          KK, LS, V, WLAM, RES
            CALL X04BAF(NOUT,REC(1))
         END IF
C
   40 CONTINUE
C
      RETURN
C
99999 FORMAT (1X,'Exit H02BBF - Optimum IP solution found. ',/)
99998 FORMAT (1X,'Exit H02BBF - First IP solution found. ',/)
99997 FORMAT (1X,'Exit H02BBF - Best IP solution found. ',/)
99996 FORMAT (1X,'Final IP objective value =',G16.7)
99995 FORMAT (//1X,'Varbl',1X,'State',5X,'Value',5X,'Lower Bound',3X,
     *       'Upper Bound',4X,'Lagr Mult',3X,'Residual',/)
99994 FORMAT (//1X,'L Con',1X,'State',5X,'Value',5X,'Lower Bound',3X,
     *       'Upper Bound',4X,'Lagr Mult',3X,'Residual',/)
99993 FORMAT (1X,'V',I3,4X,A2,1X,1P,3G14.6,1P,2G12.4)
99992 FORMAT (1X,'V',I3,4X,A2,1X,1P,G14.6,'     None     ',G14.6,1P,
     *       2G12.4)
99991 FORMAT (1X,'V',I3,4X,A2,1X,1P,2G14.6,'     None     ',1P,2G12.4)
99990 FORMAT (1X,'V',I3,4X,A2,1X,1P,G14.6,
     *       '     None          None     ',1P,2G12.4)
99989 FORMAT (1X,'L',I3,4X,A2,1X,1P,G14.6,'     None     ',G14.6,1P,
     *       2G12.4)
99988 FORMAT (1X,'L',I3,4X,A2,1X,1P,2G14.6,'     None     ',1P,2G12.4)
99987 FORMAT (1X,'L',I3,4X,A2,1X,1P,G14.6,
     *       '     None          None     ',1P,2G12.4)
99986 FORMAT (1X,'L',I3,4X,A2,1X,1P,3G14.6,1P,2G12.4)
99985 FORMAT (//1X,'Total of ',I5,' nodes investigated.',/)
      END
