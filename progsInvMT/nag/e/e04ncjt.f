      SUBROUTINE E04NCJ(PRBTYP,ISDEL,ITER,JADD,JDEL,MSGLVL,NACTIV,NFREE,
     *                  N,NCLIN,NRANK,LDR,LDT,NZ,NRZ,ISTATE,ALFA,CONDRZ,
     *                  CONDT,GZRNRM,NUMINF,SUMINF,CTX,SSQ,AX,R,T,X,
     *                  WORK)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1062 (JUL 1993).
C     MARK 17 REVISED. IER-1574 (JUN 1995).
C
C     ******************************************************************
C     E04NCJ  prints various levels of output for  E04NCZ.
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
C        ge  20        constraint status,  x  and  Ax.
C
C        ge  30        diagonals of  T  and  R.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October-1984.
C     This version of E04NCJ dated 21-Oct-1992.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           MLINE1, MLINE2
      PARAMETER         (MLINE1=50000,MLINE2=50000)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, CONDRZ, CONDT, CTX, GZRNRM, SSQ, SUMINF
      INTEGER           ISDEL, ITER, JADD, JDEL, LDR, LDT, MSGLVL, N,
     *                  NACTIV, NCLIN, NFREE, NRANK, NRZ, NUMINF, NZ
      CHARACTER*2       PRBTYP
C     .. Array Arguments ..
      DOUBLE PRECISION  AX(*), R(LDR,*), T(LDT,*), WORK(N), X(N)
      INTEGER           ISTATE(*)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  OBJ
      INTEGER           I, ITN, J, K, KADD, KDEL, NART, NDF
      LOGICAL           FIRST, LINOBJ, NEWSET, PRTHDR
      CHARACTER*2       LADD, LDEL
C     .. Local Arrays ..
      CHARACTER*2       LSTATE(0:5)
      CHARACTER*120     REC(4)
C     .. External Subroutines ..
      EXTERNAL          DCOPY, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MIN, MOD
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
C     .. Data statements ..
      DATA              LSTATE(0), LSTATE(1), LSTATE(2)/'  ', 'L ',
     *                  'U '/
      DATA              LSTATE(3), LSTATE(4), LSTATE(5)/'E ', 'F ',
     *                  'A '/
C     .. Executable Statements ..
      IF (MSGLVL.GE.15) THEN
         IF (ISUMM.GE.0) THEN
            WRITE (REC,FMT=99999) PRBTYP, ITER
            CALL X04BAY(ISUMM,4,REC)
         END IF
      END IF
C
      IF (MSGLVL.GE.5) THEN
C
         FIRST = ITER .EQ. 0
         LINOBJ = NRANK .EQ. 0
C
         ITN = MOD(ITER,10000)
         NDF = MOD(NRZ,10000)
C
         NART = NZ - NRZ
C
         IF (JDEL.GT.0) THEN
            KDEL = ISDEL
         ELSE IF (JDEL.LT.0) THEN
            JDEL = NART + 1
            KDEL = 5
         ELSE
            KDEL = 0
         END IF
C
         IF (JADD.GT.0) THEN
            KADD = ISTATE(JADD)
         ELSE
            KADD = 0
         END IF
C
         LDEL = LSTATE(KDEL)
         LADD = LSTATE(KADD)
C
         IF (NUMINF.GT.0) THEN
            OBJ = SUMINF
         ELSE
            OBJ = SSQ + CTX
         END IF
C        ---------------------------------------------------------------
C        If necessary, print a header.
C        Print a single line of information.
C        ---------------------------------------------------------------
         IF (ISUMM.GE.0) THEN
C           -----------------------------------
C           Terse line for the Monitoring file.
C           -----------------------------------
            NEWSET = LINES1 .GE. MLINE1
            PRTHDR = MSGLVL .GE. 15 .OR. FIRST .OR. NEWSET
C
            IF (PRTHDR) THEN
               IF (LINOBJ) THEN
                  WRITE (REC,FMT=99998)
                  CALL X04BAY(ISUMM,3,REC)
               ELSE
                  WRITE (REC,FMT=99996)
                  CALL X04BAY(ISUMM,3,REC)
               END IF
               LINES1 = 0
            END IF
C
            IF (LINOBJ) THEN
               WRITE (REC,FMT=99994) ITN, JDEL, LDEL, JADD, LADD, ALFA,
     *           NUMINF, OBJ, N - NFREE, NACTIV, NART, NDF, GZRNRM,
     *           CONDT
               CALL X04BAF(ISUMM,REC(1))
            ELSE
               WRITE (REC,FMT=99994) ITN, JDEL, LDEL, JADD, LADD, ALFA,
     *           NUMINF, OBJ, N - NFREE, NACTIV, NART, NDF, GZRNRM,
     *           CONDT, CONDRZ
               CALL X04BAF(ISUMM,REC(1))
            END IF
            LINES1 = LINES1 + 1
         END IF
C
         IF (IPRINT.GE.0 .AND. ISUMM.NE.IPRINT) THEN
C           ------------------------------
C           Terse line for the Print file.
C           ------------------------------
            NEWSET = LINES2 .GE. MLINE2
            PRTHDR = FIRST .OR. NEWSET
C
            IF (PRTHDR) THEN
               WRITE (REC,FMT=99997)
               CALL X04BAY(IPRINT,3,REC)
               LINES2 = 0
            END IF
C
            WRITE (REC,FMT=99995) ITN, ALFA, NUMINF, OBJ, GZRNRM
            CALL X04BAF(IPRINT,REC(1))
            LINES2 = LINES2 + 1
         END IF
C
         IF (MSGLVL.GE.20) THEN
            IF (ISUMM.GE.0) THEN
               WRITE (REC,FMT=99993) PRBTYP
               CALL X04BAY(ISUMM,3,REC)
               WRITE (REC,FMT=99992)
               CALL X04BAY(ISUMM,2,REC)
               DO 20 I = 1, N, 5
                  WRITE (REC,FMT=99987) (X(J),ISTATE(J),J=I,MIN(I+4,N))
                  CALL X04BAF(ISUMM,REC(1))
   20          CONTINUE
               IF (NCLIN.GT.0) THEN
                  WRITE (REC,FMT=99991)
                  CALL X04BAY(ISUMM,2,REC)
                  DO 40 I = 1, NCLIN, 5
                     WRITE (REC,FMT=99987) (AX(K),ISTATE(N+K),K=I,
     *                 MIN(I+4,NCLIN))
                     CALL X04BAF(ISUMM,REC(1))
   40             CONTINUE
               END IF
C
               IF (MSGLVL.GE.30) THEN
C                 ------------------------------------------------------
C                 Print the diagonals of  T  and  R.
C                 ------------------------------------------------------
                  IF (NACTIV.GT.0) THEN
                     CALL DCOPY(NACTIV,T(NACTIV,NZ+1),LDT-1,WORK,1)
                     WRITE (REC,FMT=99990) PRBTYP
                     CALL X04BAY(ISUMM,2,REC)
                     DO 60 I = 1, NACTIV, 5
                        WRITE (REC,FMT=99986) (WORK(J),J=I,
     *                    MIN(I+4,NACTIV))
                        CALL X04BAF(ISUMM,REC(1))
   60                CONTINUE
                  END IF
                  IF (NRANK.GT.0) THEN
                     WRITE (REC,FMT=99989) PRBTYP
                     CALL X04BAY(ISUMM,2,REC)
                     DO 80 I = 1, NRANK, 5
                        WRITE (REC,FMT=99986) (R(J,J),J=I,MIN(I+4,NRANK)
     *                    )
                        CALL X04BAF(ISUMM,REC(1))
   80                CONTINUE
                  END IF
               END IF
               WRITE (REC,FMT=99988)
               CALL X04BAY(ISUMM,3,REC)
            END IF
         END IF
      END IF
C
      RETURN
C
C
C     End of  E04NCJ. (LSPRT)
C
99999 FORMAT (//' ',A2,' iteration',I5,/' =================')
99998 FORMAT (//' Itn Jdel  Jadd      Step Ninf  Sinf/Objective  Bnd  ',
     *       'Lin  Art   Zr  Norm Gz   Cond T')
99997 FORMAT (//' Itn     Step Ninf Sinf/Objective  Norm Gz')
99996 FORMAT (//' Itn Jdel  Jadd      Step Ninf  Sinf/Objective  Bnd  ',
     *       'Lin  Art   Zr  Norm Gz   Cond T  Cond Rz')
99995 FORMAT (I4,1P,D9.1,I5,D15.6,D9.1)
99994 FORMAT (I4,I5,A1,I5,A1,1P,D9.1,I5,D16.8,4I5,3D9.1)
99993 FORMAT (/' Values and status of the ',A2,' constraints',/' -----',
     *       '----------------------------------')
99992 FORMAT (/' Variables...')
99991 FORMAT (/' General linear constraints...')
99990 FORMAT (/' Diagonals of ',A2,' working set factor T')
99989 FORMAT (/' Diagonals of ',A2,' triangle R         ')
99988 FORMAT (//' ----------------------------------------------------',
     *       '-------------------------------------------')
99987 FORMAT (1X,5(1P,D15.6,I5))
99986 FORMAT (1P,5D15.6)
      END
