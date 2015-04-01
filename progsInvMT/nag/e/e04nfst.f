      SUBROUTINE E04NFS(PRBTYP,HEADER,RSET,MSGLVL,ITER,ISDEL,JDEL,JADD,
     *                  N,NCLIN,NACTIV,NFREE,NZ,NRZ,LDR,LDT,ISTATE,ALFA,
     *                  CONDRZ,CONDT,DRZZ,GRZNRM,NUMINF,SUMINF,NOTOPT,
     *                  OBJQP,TRUSML,AX,R,T,X,WORK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1592 (JUN 1995).
C
C     ==================================================================
C     E04NFS prints various levels of output for E04NFZ.
C
C           msg        cumulative result
C           ---        -----------------
C
C       .le.  0        no output.
C
C       .eq.  1        nothing now (but full output later).
C
C       .eq.  5        one terse line of output.
C
C       .ge. 10        same as 5 (but full output later).
C
C       .ge. 20        constraint status,  x  and  Ax.
C
C       .ge. 30        diagonals of  T  and  R.
C
C
C     Original version of E04NFS written by PEG, 31-October-1984.
C     This version of  E04NFS  dated  11-Nov-92.
C     ==================================================================
C
C     .. Parameters ..
      INTEGER           MLINE1, MLINE2
      PARAMETER         (MLINE1=50000,MLINE2=50000)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, CONDRZ, CONDT, DRZZ, GRZNRM, OBJQP,
     *                  SUMINF, TRUSML
      INTEGER           ISDEL, ITER, JADD, JDEL, LDR, LDT, MSGLVL, N,
     *                  NACTIV, NCLIN, NFREE, NOTOPT, NRZ, NUMINF, NZ
      LOGICAL           HEADER, RSET
      CHARACTER*2       PRBTYP
C     .. Array Arguments ..
      DOUBLE PRECISION  AX(*), R(LDR,*), T(LDT,*), WORK(N), X(N)
      INTEGER           ISTATE(*)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  OBJ
      INTEGER           ITN, J, JJ, K, KADD, KDEL, KK, NDF
      LOGICAL           NEWSET, PRTHDR
      CHARACTER*2       LADD, LDEL
      CHARACTER*9       LCONDR, LRZZ
      CHARACTER*15      LMCHAR
C     .. Local Arrays ..
      CHARACTER*2       LSTATE(0:5)
      CHARACTER*120     REC(5)
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
C
      IF (MSGLVL.GE.15) THEN
         IF (ISUMM.GE.0) THEN
            WRITE (REC,FMT=99999) PRBTYP, ITER
            CALL X04BAY(ISUMM,5,REC)
         END IF
      END IF
C
      IF (MSGLVL.GE.5) THEN
C        ---------------------------------------------------------------
C        Some printing required.  Set up information for the terse line.
C        ---------------------------------------------------------------
         ITN = MOD(ITER,10000)
         NDF = MOD(NRZ,10000)
C
         IF (JDEL.NE.0) THEN
            IF (NOTOPT.GT.0) THEN
               WRITE (LMCHAR,FMT='( I5, 1P,D10.2 )') NOTOPT, TRUSML
            ELSE
               WRITE (LMCHAR,FMT='( 5X, 1P,D10.2 )') TRUSML
            END IF
C
            IF (JDEL.GT.0) THEN
               KDEL = ISDEL
C
            ELSE IF (JDEL.LT.0) THEN
               JDEL = NZ - NRZ + 1
               KDEL = 5
            END IF
         ELSE
            JDEL = 0
            KDEL = 0
            LMCHAR = '               '
         END IF
C
         LRZZ = '         '
         LCONDR = '         '
         IF (RSET .AND. NRZ.GT.0) THEN
            WRITE (LCONDR,FMT='( 1P,D9.1 )') CONDRZ
            IF (DRZZ.NE.ONE) WRITE (LRZZ,FMT='( 1P,D9.1 )') DRZZ
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
            OBJ = OBJQP
         END IF
C
C        ---------------------------------------------------------------
C        If necessary, print a header.
C        Print a single line of information.
C        ---------------------------------------------------------------
         IF (ISUMM.GE.0) THEN
C           -----------------------------------
C           Terse line for the Monitoring file.
C           -----------------------------------
            NEWSET = LINES1 .GE. MLINE1
            PRTHDR = MSGLVL .GE. 15 .OR. HEADER .OR. NEWSET
C
            IF (PRTHDR) THEN
               IF (PRBTYP.EQ.'QP') THEN
                  WRITE (REC,FMT=99997)
                  CALL X04BAY(ISUMM,3,REC)
               ELSE
                  WRITE (REC,FMT=99996)
                  CALL X04BAY(ISUMM,3,REC)
               END IF
               LINES1 = 0
            END IF
C
            WRITE (REC,FMT=99994) ITN, JDEL, LDEL, JADD, LADD, ALFA,
     *        NUMINF, OBJ, N - NFREE, NACTIV, NZ - NRZ, NDF, GRZNRM,
     *        LMCHAR, CONDT, LCONDR, LRZZ
            CALL X04BAF(ISUMM,REC(1))
            LINES1 = LINES1 + 1
         END IF
C
         IF (IPRINT.GE.0 .AND. ISUMM.NE.IPRINT) THEN
C           ------------------------------
C           Terse line for the Print file.
C           ------------------------------
            NEWSET = LINES2 .GE. MLINE2
            PRTHDR = HEADER .OR. NEWSET
C
            IF (PRTHDR) THEN
               WRITE (REC,FMT=99998)
               CALL X04BAY(IPRINT,3,REC)
               LINES2 = 0
            END IF
C
            WRITE (REC,FMT=99995) ITN, ALFA, NUMINF, OBJ, GRZNRM
            CALL X04BAF(IPRINT,REC(1))
C
            LINES2 = LINES2 + 1
         END IF
C
         IF (MSGLVL.GE.20) THEN
            IF (ISUMM.GE.0) THEN
               WRITE (REC,FMT=99993) PRBTYP
               CALL X04BAY(ISUMM,3,REC)
               WRITE (REC,FMT=99992)
               CALL X04BAY(ISUMM,2,REC)
               DO 20 J = 1, N, 5
                  WRITE (REC,FMT=99991) (X(JJ),ISTATE(JJ),JJ=J,
     *              MIN(J+4,N))
                  CALL X04BAF(ISUMM,REC(1))
   20          CONTINUE
               IF (NCLIN.GT.0) THEN
                  WRITE (REC,FMT=99990)
                  CALL X04BAY(ISUMM,2,REC)
                  DO 40 K = 1, NCLIN, 5
                     WRITE (REC,FMT=99991) (AX(KK),ISTATE(N+KK),KK=K,
     *                 MIN(K+4,NCLIN))
                     CALL X04BAF(ISUMM,REC(1))
   40             CONTINUE
               END IF
C
               IF (MSGLVL.GE.30) THEN
C                 ------------------------------------------------------
C                 Print the diagonals of  T  and  R.
C                 ------------------------------------------------------
                  IF (NACTIV.GT.0) THEN
                     CALL DCOPY(NACTIV,T(1,NZ+1),LDT+1,WORK,1)
                     WRITE (REC,FMT=99989) PRBTYP
                     CALL X04BAY(ISUMM,2,REC)
                     DO 60 J = 1, NACTIV, 5
                        WRITE (REC,FMT=99988) (WORK(JJ),JJ=J,
     *                    MIN(J+4,NACTIV))
                        CALL X04BAF(ISUMM,REC(1))
   60                CONTINUE
                  END IF
                  IF (RSET .AND. NRZ.GT.0) THEN
                     WRITE (REC,FMT=99987) PRBTYP
                     CALL X04BAY(ISUMM,2,REC)
                     DO 80 J = 1, NRZ, 5
                        WRITE (REC,FMT=99988) (R(JJ,JJ),JJ=J,
     *                    MIN(J+4,NRZ))
                        CALL X04BAF(ISUMM,REC(1))
   80                CONTINUE
                  END IF
               END IF
               WRITE (REC,FMT=99986)
               CALL X04BAY(ISUMM,4,REC)
            END IF
         END IF
      END IF
C
      HEADER = .FALSE.
      JDEL = 0
      JADD = 0
      ALFA = ZERO
C
      RETURN
C
C
C     End of E04NFS.  (QPPRNT)
C
99999 FORMAT (///' ',A2,' iteration',I5,/' =================')
99998 FORMAT (//'  Itn     Step Ninf Sinf/Objective  Norm Gz')
99997 FORMAT (//'  Itn Jdel  Jadd      Step Ninf  Sinf/Objective  Bnd ',
     *       ' Lin  Art   Zr  Norm Gz NOpt    Min Lm   Cond T  Cond Rz',
     *       '     Rzz')
99996 FORMAT (//'  Itn Jdel  Jadd      Step Ninf  Sinf/Objective  Bnd ',
     *       ' Lin  Art   Zr  Norm Gz NOpt    Min Lm   Cond T')
99995 FORMAT (I5,1P,D9.1,I5,D15.6,D9.1)
99994 FORMAT (I5,I5,A1,I5,A1,1P,D9.1,I5,D16.8,4I5,D9.1,A15,D9.1,2A9)
99993 FORMAT (/' Values and status of the ',A2,' constraints',/' -----',
     *       '----------------------------------')
99992 FORMAT (/' Variables...')
99991 FORMAT (1X,5(1P,D15.6,I5))
99990 FORMAT (/' General linear constraints...')
99989 FORMAT (/' Diagonals of ',A2,' working set factor T')
99988 FORMAT (1P,5D15.6)
99987 FORMAT (/' Diagonals of ',A2,' triangle Rz        ')
99986 FORMAT (///' ---------------------------------------------------',
     *       '--------------------------------------------')
      END
