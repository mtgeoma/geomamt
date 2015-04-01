      SUBROUTINE E04UCT(KTCOND,CONVRG,LSUMRY,MSGNP,MSGQP,LDR,LDT,N,
     *                  NCLIN,NCNLN,NCTOTL,NACTIV,LINACT,NLNACT,NZ,
     *                  NFREE,MAJIT0,MAJITS,MINITS,ISTATE,ALFA,NFUN,
     *                  CONDHZ,CONDH,CONDT,OBJALF,OBJF,GZNORM,CVNORM,AX,
     *                  C,R,T,VIOLN,X,WORK)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-723 (DEC 1989).
C     MARK 16 REVISED. IER-1088 (JUL 1993).
C     MARK 17 REVISED. IER-1608 (JUN 1995).
C
C     ==================================================================
C     E04UCT  prints various levels of output for E04UCZ and E04UPZ.
C
C           Msg        Cumulative result
C           ---        -----------------
C
C        le   0        no output.
C
C        eq   1        nothing now (but full output later).
C
C        eq   5        one terse line of output.
C
C        ge  10        same as 5 (but full output later).
C
C        ge  20        objective function,  x,  Ax  and  c.
C
C        ge  30        diagonals of  T  and  R.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original Fortran 66 version written November-1982.
C     This version of  E04UCT  dated  21-Oct-92.
C     ==================================================================
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, CONDH, CONDHZ, CONDT, CVNORM, GZNORM,
     *                  OBJALF, OBJF
      INTEGER           LDR, LDT, LINACT, MAJIT0, MAJITS, MINITS, MSGNP,
     *                  MSGQP, N, NACTIV, NCLIN, NCNLN, NCTOTL, NFREE,
     *                  NFUN, NLNACT, NZ
      LOGICAL           CONVRG
      CHARACTER*5       LSUMRY
C     .. Array Arguments ..
      DOUBLE PRECISION  AX(*), C(*), R(LDR,*), T(LDT,*), VIOLN(*),
     *                  WORK(N), X(N)
      INTEGER           ISTATE(NCTOTL)
      LOGICAL           KTCOND(2)
C     .. Scalars in Common ..
      DOUBLE PRECISION  RHODMP, RHOMAX, RHONRM, SCALE
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           INCRUN
C     .. Local Scalars ..
      DOUBLE PRECISION  CVIOLS
      INTEGER           I, INCT, J, K, MJR, MNR, NDF, NEVAL
      LOGICAL           FIRST, NEWSET, NLNCON, PRTHDR
C     .. Local Arrays ..
      CHARACTER*132     REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      EXTERNAL          DNRM2
C     .. External Subroutines ..
      EXTERNAL          DCOPY, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MIN, MOD
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /DE04UC/RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
C     .. Executable Statements ..
C
      IF (MSGNP.GE.20) THEN
         IF (ISUMM.GE.0) THEN
            WRITE (REC,FMT=99999) MAJITS
            CALL X04BAY(ISUMM,4,REC)
         END IF
      END IF
C
      IF (MSGNP.GE.5) THEN
C
         MJR = MOD(MAJITS,1000)
         MNR = MOD(MINITS,1000)
         NEVAL = MOD(NFUN,1000)
         NDF = MOD(NZ,1000)
         NLNCON = NCNLN .GT. 0
         FIRST = MAJITS .EQ. MAJIT0
C
C        ---------------------------------------------------------------
C        If necessary, print a header.
C        Print a single line of information.
C        ---------------------------------------------------------------
         IF (ISUMM.GE.0) THEN
C           -----------------------------------
C           Terse line for the Monitoring file.
C           -----------------------------------
            NEWSET = LINES1 .GE. 50000
            PRTHDR = MSGQP .GT. 0 .OR. FIRST .OR. MSGNP .GE. 20 .OR.
     *               NEWSET
C
            IF (PRTHDR) THEN
               IF (NLNCON) THEN
                  WRITE (REC,FMT=99998)
                  CALL X04BAY(ISUMM,3,REC)
               ELSE
                  WRITE (REC,FMT=99996)
                  CALL X04BAY(ISUMM,3,REC)
               END IF
               LINES1 = 0
            END IF
C
            IF (NLNCON) THEN
               WRITE (REC,FMT=99997) MJR, MNR, ALFA, NEVAL, OBJALF,
     *           GZNORM, CVNORM, NDF, N - NFREE, LINACT, NLNACT,
     *           SCALE*RHONRM, CONDH, CONDHZ, CONDT, CONVRG, KTCOND(1),
     *           KTCOND(2), LSUMRY
               CALL X04BAF(ISUMM,REC(1))
            ELSE
               WRITE (REC,FMT=99995) MJR, MNR, ALFA, NEVAL, OBJALF,
     *           GZNORM, NDF, N - NFREE, LINACT, CONDH, CONDHZ, CONDT,
     *           CONVRG, KTCOND(1), KTCOND(2), LSUMRY
               CALL X04BAF(ISUMM,REC(1))
            END IF
            LINES1 = LINES1 + 1
         END IF
C
         IF (IPRINT.GE.0 .AND. ISUMM.NE.IPRINT) THEN
C           ------------------------------
C           Terse line for the Print file.
C           ------------------------------
            NEWSET = LINES2 .GE. 50000
            PRTHDR = MSGQP .GT. 0 .OR. FIRST .OR. NEWSET
C
            IF (PRTHDR) THEN
               IF (NLNCON) THEN
                  WRITE (REC,FMT=99994)
                  CALL X04BAY(IPRINT,3,REC)
               ELSE
                  WRITE (REC,FMT=99992)
                  CALL X04BAY(IPRINT,3,REC)
               END IF
               LINES2 = 0
            END IF
C
            IF (NLNCON) THEN
               WRITE (REC,FMT=99993) MJR, MNR, ALFA, OBJALF, GZNORM,
     *           CVNORM, CONDHZ, LSUMRY
               CALL X04BAF(IPRINT,REC(1))
            ELSE
               WRITE (REC,FMT=99991) MJR, MNR, ALFA, OBJALF, GZNORM,
     *           CONDHZ, LSUMRY
               CALL X04BAF(IPRINT,REC(1))
            END IF
            LINES2 = LINES2 + 1
         END IF
C
         IF (MSGNP.GE.20) THEN
            IF (ISUMM.GE.0) THEN
               IF (NCNLN.EQ.0) THEN
                  WRITE (REC,FMT=99990) OBJF
                  CALL X04BAY(ISUMM,2,REC)
               ELSE
                  CVIOLS = DNRM2(NCNLN,VIOLN,1)
                  WRITE (REC,FMT=99989) OBJF, CVIOLS
                  CALL X04BAY(ISUMM,2,REC)
               END IF
C
C              ---------------------------------------------------------
C              Print the constraint values.
C              ---------------------------------------------------------
               WRITE (REC,FMT=99988)
               CALL X04BAY(ISUMM,3,REC)
               WRITE (REC,FMT=99987)
               CALL X04BAY(ISUMM,2,REC)
               DO 20 I = 1, N, 5
                  WRITE (REC,FMT=99982) (X(J),ISTATE(J),J=I,MIN(I+4,N))
                  CALL X04BAF(ISUMM,REC(1))
   20          CONTINUE
               IF (NCLIN.GT.0) THEN
                  WRITE (REC,FMT=99986)
                  CALL X04BAY(ISUMM,2,REC)
                  DO 40 I = 1, NCLIN, 5
                     WRITE (REC,FMT=99982) (AX(K),ISTATE(N+K),K=I,
     *                 MIN(I+4,NCLIN))
                     CALL X04BAF(ISUMM,REC(1))
   40             CONTINUE
               END IF
               IF (NCNLN.GT.0) THEN
                  WRITE (REC,FMT=99985)
                  CALL X04BAY(ISUMM,2,REC)
                  DO 60 I = 1, NCNLN, 5
                     WRITE (REC,FMT=99982) (C(K),ISTATE(N+NCLIN+K),K=I,
     *                 MIN(I+4,NCNLN))
                     CALL X04BAF(ISUMM,REC(1))
   60             CONTINUE
               END IF
C
               IF (MSGNP.GE.30) THEN
C                 ------------------------------------------------------
C                 Print the diagonals of  T  and  R.
C                 ------------------------------------------------------
                  INCT = LDT - 1
                  IF (NACTIV.GT.0) THEN
                     CALL DCOPY(NACTIV,T(NACTIV,NZ+1),INCT,WORK,1)
                     WRITE (REC,FMT=99984)
                     CALL X04BAY(ISUMM,2,REC)
                     DO 80 I = 1, NACTIV, 5
                        WRITE (REC,FMT=99981) (WORK(J),J=I,
     *                    MIN(I+4,NACTIV))
                        CALL X04BAF(ISUMM,REC(1))
   80                CONTINUE
                  END IF
                  WRITE (REC,FMT=99983)
                  CALL X04BAY(ISUMM,2,REC)
                  DO 100 I = 1, N, 5
                     WRITE (REC,FMT=99981) (R(J,J),J=I,MIN(I+4,N))
                     CALL X04BAF(ISUMM,REC(1))
  100             CONTINUE
               END IF
            END IF
         END IF
      END IF
C
      IF (MSGNP.GE.20) THEN
         IF (ISUMM.GE.0) THEN
            WRITE (REC,FMT=99980)
            CALL X04BAY(ISUMM,3,REC)
         END IF
      END IF
C
      LSUMRY(1:2) = '  '
      LSUMRY(4:5) = '  '
C
      RETURN
C
C
C     End of E04UCT. (NPPRT)
C
99999 FORMAT (//' Major iteration',I5,/' ====================')
99998 FORMAT (//'  Maj  Mnr    Step Nfun  Merit Function Norm Gz  Viol',
     *       'tn   Nz  Bnd  Lin  Nln Penalty  Cond H Cond Hz  Cond T C',
     *       'onv')
99997 FORMAT (2I5,1P,D8.1,I5,D16.8,2D8.1,4I5,4D8.1,1X,L1,1X,2L1,A5)
99996 FORMAT (//'  Maj  Mnr    Step Nfun       Objective Norm Gz   Nz ',
     *       ' Bnd  Lin  Cond H Cond Hz  Cond T Conv')
99995 FORMAT (2I5,1P,D8.1,I5,D16.8,D8.1,3I5,3D8.1,1X,L1,1X,2L1,A5)
99994 FORMAT (//'  Maj  Mnr    Step Merit Function Norm Gz  Violtn Con',
     *       'd Hz')
99993 FORMAT (2I5,1P,D8.1,D15.6,3D8.1,2X,A5)
99992 FORMAT (//'  Maj  Mnr    Step      Objective Norm Gz Cond Hz')
99991 FORMAT (2I5,1P,D8.1,D15.6,2D8.1,2X,A5)
99990 FORMAT (/' Nonlinear objective value = ',1P,D15.6)
99989 FORMAT (/' Nonlinear objective value = ',1P,D15.6,'   Norm of th',
     *       'e nonlinear constraint violations = ',D15.6)
99988 FORMAT (/' Values of the constraints and their predicted status',
     *       /' ----------------------------------------------------')
99987 FORMAT (/' Variables                  ')
99986 FORMAT (/' General linear constraints ')
99985 FORMAT (/' Nonlinear constraints      ')
99984 FORMAT (/' Diagonals of  T  =         ')
99983 FORMAT (/' Diagonals of  R  =         ')
99982 FORMAT (1X,5(1P,D15.6,I4))
99981 FORMAT (1P,5D15.6)
99980 FORMAT (' ======================================================',
     *       '===================================',//)
      END
