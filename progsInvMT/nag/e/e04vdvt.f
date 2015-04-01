      SUBROUTINE E04VDV(NERROR,LIWORK,LWORK,LITOTL,LWTOTL,NROWA,N,NCLIN,
     *                  NCNLN,NCTOTL,ISTATE,KACTIV,LCRASH,NAMED,NAMES,
     *                  LENNAM,BIGBND,A,BL,BU,FEATOL,X)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12A REVISED. IER-507 (AUG 1986).
C
C *********************************************************************
C     E04VDV  CHECKS THE DATA INPUT TO THE VARIOUS OPTIMIZERS.
C
C     THE FOLLOWING QUANTITIES ARE NOT CHECKED...
C     NROWA, N, NCLIN, NCTOTL
C     KACTIV
C     A, X
C
C     SYSTEMS OPTIMIZATION LABORATORY, STANFORD UNIVERSITY.
C     VERSION OF APRIL 1982.  REV. OCT. 1982.
C *********************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND
      INTEGER           LCRASH, LENNAM, LITOTL, LIWORK, LWORK, LWTOTL,
     *                  N, NCLIN, NCNLN, NCTOTL, NERROR, NROWA
      LOGICAL           NAMED
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NROWA,N), BL(NCTOTL), BU(NCTOTL),
     *                  FEATOL(NCTOTL), X(N)
      INTEGER           ISTATE(NCTOTL), KACTIV(N), NAMES(4,LENNAM)
C     .. Scalars in Common ..
      INTEGER           ISTART, MSG, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  B1, B2, FTOL, ONE, TEST, ZERO
      INTEGER           IS, J, K, L, L1, L2, NPLIN
      LOGICAL           OK
C     .. Local Arrays ..
      CHARACTER*2       ID(9)
      CHARACTER*95      REC(3)
C     .. External Subroutines ..
      EXTERNAL          X04BAF
C     .. Common blocks ..
      COMMON            /AE04VC/NOUT, MSG, ISTART
C     .. Data statements ..
      DATA              ID(1), ID(2), ID(3), ID(4), ID(5)/'VA', 'RB',
     *                  'L ', 'LN', 'CO'/
      DATA              ID(6), ID(7), ID(8), ID(9)/'N ', 'NL', 'CO',
     *                  'N '/
      DATA              ONE/1.0D+0/, ZERO/0.0D+0/
C     .. Executable Statements ..
C
      NERROR = 0
C
C ---------------------------------------------------------------------
C     CHECK THAT THERE IS ENOUGH WORKSPACE TO SOLVE THE PROBLEM.
C ---------------------------------------------------------------------
      OK = LITOTL .LE. LIWORK .AND. LWTOTL .LE. LWORK
      IF (OK .AND. MSG.LE.0) GO TO 20
      IF (MSG.GE.0) THEN
         WRITE (REC,FMT=99999) LIWORK, LWORK, LITOTL, LWTOTL
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         CALL X04BAF(NOUT,REC(3))
      END IF
      IF (OK) GO TO 20
      NERROR = NERROR + 1
      IF (MSG.GE.0) THEN
         WRITE (REC,FMT=99998)
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
      END IF
C
C ---------------------------------------------------------------------
C     CHECK THE BOUNDS ON ALL VARIABLES AND CONSTRAINTS.
C ---------------------------------------------------------------------
   20 DO 40 J = 1, NCTOTL
         B1 = BL(J)
         B2 = BU(J)
         OK = B1 .LE. B2
         IF (OK) GO TO 40
         NERROR = NERROR + 1
         K = J
         L1 = 1
         IF (J.GT.N) K = J - N
         IF (J.GT.N) L1 = 4
         IF (J.GT.N+NCLIN) K = K - NCLIN
         IF (J.GT.N+NCLIN) L1 = 7
         L2 = L1 + 2
         IF ( .NOT. NAMED .AND. MSG.GE.0) THEN
            WRITE (REC,FMT=99997) (ID(L),L=L1,L2), K, B1, B2
            CALL X04BAF(NOUT,REC(1))
            CALL X04BAF(NOUT,REC(2))
         END IF
         IF (NAMED .AND. MSG.GE.0) THEN
            WRITE (REC,FMT=99996) (NAMES(L,J),L=1,4), B1, B2
            CALL X04BAF(NOUT,REC(1))
            CALL X04BAF(NOUT,REC(2))
         END IF
   40 CONTINUE
C
C ---------------------------------------------------------------------
C     CHECK  BIGBND  AND  FEATOL.
C ---------------------------------------------------------------------
      OK = BIGBND .GT. ZERO
      IF (OK) GO TO 60
      NERROR = NERROR + 1
      IF (MSG.GE.0) THEN
         WRITE (REC,FMT=99995) BIGBND
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
      END IF
C
   60 DO 80 J = 1, NCTOTL
         FTOL = FEATOL(J)
         TEST = ONE + FTOL
         OK = TEST .GT. ONE
         IF (OK) GO TO 80
         IF (MSG.GE.0) THEN
            WRITE (REC,FMT=99994) J, FTOL
            CALL X04BAF(NOUT,REC(1))
            CALL X04BAF(NOUT,REC(2))
         END IF
   80 CONTINUE
C
C ---------------------------------------------------------------------
C     IF WARM START, CHECK  ISTATE.
C ---------------------------------------------------------------------
  100 IF (LCRASH.EQ.0) GO TO 140
      NPLIN = N + NCLIN
C
      DO 120 J = 1, NPLIN
         IS = ISTATE(J)
         OK = IS .GE. (-2) .AND. IS .LE. 4
         IF (OK) GO TO 120
         NERROR = NERROR + 1
         IF (MSG.GE.0) THEN
            WRITE (REC,FMT=99993) J, IS
            CALL X04BAF(NOUT,REC(1))
            CALL X04BAF(NOUT,REC(2))
         END IF
  120 CONTINUE
C
  140 RETURN
C
C
C     END OF E04VDV  ( CHKDAT )
99999 FORMAT (/' WORKSPACE PROVIDED IS     IW(',I6,'),  W(',I6,').',
     *  /' TO SOLVE PROBLEM WE NEED  IW(',I6,'),  W(',I6,').')
99998 FORMAT (/' XXX  NOT ENOUGH WORKSPACE TO SOLVE PROBLEM.')
99997 FORMAT (/' XXX  THE BOUNDS ON  ',2A2,A1,I3,'  ARE INCONSISTENT. ',
     *  '  BL =',G16.7,'   BU =',G16.7)
99996 FORMAT (/' XXX  THE BOUNDS ON  ',4A2,'  ARE INCONSISTENT.   BL =',
     *  G16.7,'   BU =',G16.7)
99995 FORMAT (/' XXX  BIGBND  IS NOT POSITIVE...',G16.6)
99994 FORMAT (/' ***  WARNING -- FEATOL(',I4,' )  IS LESS THAN MACHINE',
     *  ' PRECISION...',G16.6)
99993 FORMAT (/' XXX  COMPONENT',I5,'  OF  ISTATE  IS OUT OF RANGE...',
     *  I10)
      END
