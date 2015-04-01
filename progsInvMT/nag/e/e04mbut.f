      SUBROUTINE E04MBU(LP,NROWA,NROWRT,NCOLRT,N,NCLIN,NCLIN0,NCTOTL,
     *                  NFREE,ISDEL,NACTIV,NCOLZ,ITER,JADD,JDEL,ALFA,
     *                  CONDT,NUMINF,SUMINF,OBJLP,ISTATE,KFREE,A,RT,X,
     *                  WRK1,WRK2)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12B REVISED. IER-536 (FEB 1987).
C
C *********************************************************************
C
C     E04MBU  PRINTS VARIOUS LEVELS OF OUTPUT FOR  E04MBY.
C
C           MSG    CUMULATIVE RESULT
C           ---    -----------------
C
C        LE   0    NO OUTPUT.
C
C        EQ   1    NOTHING NOW (BUT FULL OUTPUT LATER).
C
C        EQ   5    ONE TERSE LINE OF OUTPUT.
C
C        GE  10    SAME AS 5 (BUT FULL OUTPUT LATER).
C
C        GE  15    NOTHING MORE IF  ITER .LT. ISTART.
C                  OTHERWISE,  X,  ISTATE  AND  KFREE.
C
C        GE  20    MULTIPLIERS (PRINTED OUTSIDE E04MBU).
C                  THE ARRAY  AX.
C
C        GE  30    DIAGONALS OF  T.
C
C        GE  80    DEBUG OUTPUT.
C
C        EQ  99    A,  BL,  BU,  CVEC,  X  (CALLED FROM E04MBW).
C
C     SYSTEMS OPTIMIZATION LABORATORY, STANFORD UNIVERSITY.
C     VERSION OF DECEMBER 1981.  REV. NOV. 1982.
C *********************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, CONDT, OBJLP, SUMINF
      INTEGER           ISDEL, ITER, JADD, JDEL, N, NACTIV, NCLIN,
     *                  NCLIN0, NCOLRT, NCOLZ, NCTOTL, NFREE, NROWA,
     *                  NROWRT, NUMINF
      LOGICAL           LP
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NROWA,N), RT(NROWRT,NCOLRT), WRK1(N),
     *                  WRK2(NCLIN0), X(N)
      INTEGER           ISTATE(NCTOTL), KFREE(N)
C     .. Scalars in Common ..
      INTEGER           ISTART, MSG, NOUT
C     .. Local Scalars ..
      INTEGER           I, INCT, J, K, L1, L2, LADDI, LDELI, LENT, LROWA
      CHARACTER*1       LADD, LDEL
C     .. Local Arrays ..
      CHARACTER*1       LSTATE(5)
      CHARACTER*80      REC(6)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DCOPY, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Common blocks ..
      COMMON            /AE04VC/NOUT, MSG, ISTART
C     .. Data statements ..
      DATA              LSTATE(1), LSTATE(2)/' ', 'L'/
      DATA              LSTATE(3), LSTATE(4)/'U', 'E'/
      DATA              LSTATE(5)/'T'/
C     .. Executable Statements ..
C
      IF (MSG.LT.5) GO TO 220
C
      LDELI = 0
      LADDI = 0
      IF (JDEL.GT.0) LDELI = ISDEL
      IF (JADD.GT.0) LADDI = ISTATE(JADD)
      LDEL = LSTATE(LDELI+1)
      LADD = LSTATE(LADDI+1)
      IF (MSG.GE.15) GO TO 20
C
C ---------------------------------------------------------------------
C     PRINT HEADING (POSSIBLY) AND TERSE LINE.
C ---------------------------------------------------------------------
      IF ( .NOT. LP .AND. ITER.EQ.0) THEN
         WRITE (REC,FMT=99998)
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         CALL X04BAF(NOUT,REC(3))
      END IF
      IF (LP .AND. ITER.EQ.0) THEN
         WRITE (REC,FMT=99997)
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         CALL X04BAF(NOUT,REC(3))
      END IF
      IF ( .NOT. LP) THEN
         WRITE (REC,FMT=99996) ITER, JDEL, LDEL, JADD, LADD, ALFA,
     *     CONDT, NUMINF, SUMINF
         CALL X04BAF(NOUT,REC(1))
      END IF
      IF (LP) THEN
         WRITE (REC,FMT=99996) ITER, JDEL, LDEL, JADD, LADD, ALFA,
     *     CONDT, NUMINF, SUMINF, OBJLP
         CALL X04BAF(NOUT,REC(1))
      END IF
      GO TO 220
C
C ---------------------------------------------------------------------
C     PRINT TERSE LINE,  X,  ISTATE  AND  KFREE.
C ---------------------------------------------------------------------
   20 WRITE (REC,FMT=99999) ITER
      DO 40 I = 1, 6
         CALL X04BAF(NOUT,REC(I))
   40 CONTINUE
      IF ( .NOT. LP) THEN
         WRITE (REC,FMT=99998)
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         CALL X04BAF(NOUT,REC(3))
      END IF
      IF (LP) THEN
         WRITE (REC,FMT=99997)
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         CALL X04BAF(NOUT,REC(3))
      END IF
      IF ( .NOT. LP) THEN
         WRITE (REC,FMT=99996) ITER, JDEL, LDEL, JADD, LADD, ALFA,
     *     CONDT, NUMINF, SUMINF
         CALL X04BAF(NOUT,REC(1))
      END IF
      IF (LP) THEN
         WRITE (REC,FMT=99996) ITER, JDEL, LDEL, JADD, LADD, ALFA,
     *     CONDT, NUMINF, SUMINF, OBJLP
         CALL X04BAF(NOUT,REC(1))
      END IF
      WRITE (REC,FMT=99995)
      CALL X04BAF(NOUT,REC(1))
      CALL X04BAF(NOUT,REC(2))
      DO 60 I = 1, N, 5
         WRITE (REC,FMT=99994) (X(J),J=I,MIN(N,I+4))
         CALL X04BAF(NOUT,REC(1))
   60 CONTINUE
      WRITE (REC,FMT=99993)
      CALL X04BAF(NOUT,REC(1))
      CALL X04BAF(NOUT,REC(2))
      DO 80 I = 1, N, 10
         WRITE (REC,FMT=99992) (ISTATE(J),J=I,MIN(N,I+9))
         CALL X04BAF(NOUT,REC(1))
   80 CONTINUE
      L1 = N + 1
      L2 = N + NCLIN
      IF (L1.LE.L2) THEN
         WRITE (REC,FMT=99991)
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         DO 100 I = L1, L2, 10
            WRITE (REC,FMT=99992) (ISTATE(J),J=I,MIN(L2,I+9))
            CALL X04BAF(NOUT,REC(1))
  100    CONTINUE
      END IF
      IF (NFREE.GT.0) THEN
         WRITE (REC,FMT=99990)
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         DO 120 I = 1, NFREE, 10
            WRITE (REC,FMT=99992) (KFREE(K),K=I,MIN(NFREE,I+9))
            CALL X04BAF(NOUT,REC(1))
  120    CONTINUE
      END IF
C
C ---------------------------------------------------------------------
C     COMPUTE AND PRINT  AX.  USE  WORK = AP  TO AVOID SIDE EFFECTS.
C ---------------------------------------------------------------------
      IF (MSG.LT.20) GO TO 220
      IF (NCLIN.EQ.0) GO TO 180
      LROWA = NROWA*(N-1) + 1
      DO 140 K = 1, NCLIN
         WRK2(K) = DDOT(N,A(K,1),NROWA,X,1)
  140 CONTINUE
      WRITE (REC,FMT=99989)
      CALL X04BAF(NOUT,REC(1))
      CALL X04BAF(NOUT,REC(2))
      DO 160 I = 1, NCLIN, 5
         WRITE (REC,FMT=99994) (WRK2(K),K=I,MIN(NCLIN,I+4))
         CALL X04BAF(NOUT,REC(1))
  160 CONTINUE
C
C ---------------------------------------------------------------------
C     PRINT THE DIAGONALS OF  T.
C ---------------------------------------------------------------------
  180 IF (MSG.LT.30) GO TO 220
      LENT = NROWRT*(NACTIV-1) + 1
      INCT = NROWRT - 1
      IF (NACTIV.NE.0) CALL DCOPY(NACTIV,RT(NACTIV,NCOLZ+1)
     *                             ,INCT,WRK1,1)
      IF (NACTIV.NE.0) THEN
         WRITE (REC,FMT=99988)
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         DO 200 I = 1, NACTIV, 5
            WRITE (REC,FMT=99994) (WRK1(J),J=I,MIN(NACTIV,I+4))
            CALL X04BAF(NOUT,REC(1))
  200    CONTINUE
      END IF
C
  220 RETURN
C
C
C     END OF E04MBU  ( LPPRT )
99999 FORMAT (///' =================',/' LP ITERATION',I5,/' =========',
     *  '========')
99998 FORMAT (//'  ITN JDEL  JADD ',6X,'STEP    COND T NUMINF',8X,' SU',
     *  'MINF')
99997 FORMAT (//'  ITN JDEL  JADD ',6X,'STEP    COND T NUMINF',8X,' SU',
     *  'MINF',9X,' LPOBJ')
99996 FORMAT (I5,I5,A1,I5,A1,1P,2D10.2,I7,1P,2D15.6)
99995 FORMAT (/' LP VARIABLES')
99994 FORMAT (1P,5D15.6)
99993 FORMAT (/' STATUS OF THE LP BOUND   CONSTRAINTS')
99992 FORMAT (1X,10I4)
99991 FORMAT (/' STATUS OF THE LP GENERAL CONSTRAINTS')
99990 FORMAT (/' LIST OF FREE LP VARIABLES')
99989 FORMAT (/' VALUES OF LP GENERAL LINEAR CONSTRAINTS')
99988 FORMAT (/' DIAGONALS OF LP WORKING SET FACTOR  T  ')
      END
