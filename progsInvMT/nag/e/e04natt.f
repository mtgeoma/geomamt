      SUBROUTINE E04NAT(ORTHOG,ISDEL,ITER,JADD,JDEL,NACTIV,NCOLR,NCOLZ,
     *                  NFREE,N,NCLIN,NCLIN0,NCTOTL,NROWA,NROWRT,NCOLRT,
     *                  NHESS,ISTATE,KFREE,ALFA,CONDH,CONDT,OBJ,GFNORM,
     *                  ZTGNRM,EMAX,A,RT,X,WRK1,WRK2)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12B REVISED. IER-546 (FEB 1987).
C
C *********************************************************************
C     E04NAT  PRINTS VARIOUS LEVELS OF OUTPUT FOR  E04NAX.
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
C        GE  20    MULTIPLIERS (PRINTED OUTSIDE E04NAT).
C                  THE ARRAY  AX.
C
C        GE  30    DIAGONALS OF  T  AND  R.
C
C        GE  80    DEBUG OUTPUT.
C
C        EQ  99    CVEC  AND  HESS  (CALLED FROM E04NAV).
C
C     SYSTEMS OPTIMIZATION LABORATORY, STANFORD UNIVERSITY.
C     VERSION OF APRIL 1982.  REV. OCT. 1982.
C *********************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, CONDH, CONDT, EMAX, GFNORM, OBJ, ZTGNRM
      INTEGER           ISDEL, ITER, JADD, JDEL, N, NACTIV, NCLIN,
     *                  NCLIN0, NCOLR, NCOLRT, NCOLZ, NCTOTL, NFREE,
     *                  NHESS, NROWA, NROWRT
      LOGICAL           ORTHOG
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
      CHARACTER*1       LSTATE(6)
      CHARACTER*105     REC(6)
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
      DATA              LSTATE(6)/'V'/
C     .. Executable Statements ..
C
      IF (MSG.LT.5) GO TO 260
C
      LDELI = 0
      LADDI = 0
      IF (JDEL.GT.0) LDELI = ISDEL
      IF (JDEL.LT.0) LDELI = 5
      IF (JDEL.LT.0) JDEL = -JDEL
      IF (JADD.GT.0) LADDI = ISTATE(JADD)
      LDEL = LSTATE(LDELI+1)
      LADD = LSTATE(LADDI+1)
      IF (MSG.GE.15) GO TO 40
C
C ---------------------------------------------------------------------
C     PRINT HEADING (POSSIBLY) AND TERSE LINE.
C ---------------------------------------------------------------------
      IF (ITER.GT.0 .OR. JDEL.GT.0) GO TO 20
      IF (ORTHOG) WRITE (REC,FMT=99998)
      IF ( .NOT. ORTHOG) WRITE (REC,FMT=99997)
      CALL X04BAF(NOUT,REC(1))
      CALL X04BAF(NOUT,REC(2))
      CALL X04BAF(NOUT,REC(3))
   20 WRITE (REC,FMT=99996) ITER, JDEL, LDEL, JADD, LADD, ALFA, NHESS,
     *  OBJ, NCOLZ, GFNORM, ZTGNRM, CONDT, CONDH, EMAX
      CALL X04BAF(NOUT,REC(1))
      GO TO 260
C
C ---------------------------------------------------------------------
C     PRINT TERSE LINE,  X,  ISTATE,  KFREE.
C ---------------------------------------------------------------------
   40 WRITE (REC,FMT=99999) ITER
      DO 60 I = 1, 6
         CALL X04BAF(NOUT,REC(I))
   60 CONTINUE
      IF (ORTHOG) WRITE (REC,FMT=99998)
      IF ( .NOT. ORTHOG) WRITE (REC,FMT=99997)
      CALL X04BAF(NOUT,REC(1))
      CALL X04BAF(NOUT,REC(2))
      CALL X04BAF(NOUT,REC(3))
      WRITE (REC,FMT=99996) ITER, JDEL, LDEL, JADD, LADD, ALFA, NHESS,
     *  OBJ, NCOLZ, GFNORM, ZTGNRM, CONDT, CONDH, EMAX
      CALL X04BAF(NOUT,REC(1))
      WRITE (REC,FMT=99995)
      CALL X04BAF(NOUT,REC(1))
      CALL X04BAF(NOUT,REC(2))
      DO 80 I = 1, N, 5
         WRITE (REC,FMT=99994) (X(J),J=I,MIN(I+4,N))
         CALL X04BAF(NOUT,REC(1))
   80 CONTINUE
      WRITE (REC,FMT=99993)
      CALL X04BAF(NOUT,REC(1))
      CALL X04BAF(NOUT,REC(2))
      DO 100 I = 1, N, 10
         WRITE (REC,FMT=99992) (ISTATE(J),J=I,MIN(N,I+9))
         CALL X04BAF(NOUT,REC(1))
  100 CONTINUE
      L1 = N + 1
      L2 = N + NCLIN
      IF (L1.LE.L2) THEN
         WRITE (REC,FMT=99991)
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         DO 120 I = L1, L2, 10
            WRITE (REC,FMT=99992) (ISTATE(J),J=I,MIN(L2,I+9))
            CALL X04BAF(NOUT,REC(1))
  120    CONTINUE
      END IF
      IF (NFREE.GT.0) THEN
         WRITE (REC,FMT=99990)
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         DO 140 I = 1, NFREE, 10
            WRITE (REC,FMT=99992) (KFREE(J),J=I,MIN(NFREE,I+9))
            CALL X04BAF(NOUT,REC(1))
  140    CONTINUE
      END IF
C
C ---------------------------------------------------------------------
C     COMPUTE AND PRINT  AX.  USE  WORK  TO AVOID SIDE EFFECTS.
C ---------------------------------------------------------------------
      IF (MSG.LT.20) GO TO 260
      IF (NCLIN.EQ.0) GO TO 200
      LROWA = NROWA*(N-1) + 1
      DO 160 K = 1, NCLIN
         WRK2(K) = DDOT(N,A(K,1),NROWA,X,1)
  160 CONTINUE
      WRITE (REC,FMT=99989)
      CALL X04BAF(NOUT,REC(1))
      CALL X04BAF(NOUT,REC(2))
      DO 180 I = 1, NCLIN, 5
         WRITE (REC,FMT=99994) (WRK2(J),J=I,MIN(NCLIN,I+4))
         CALL X04BAF(NOUT,REC(1))
  180 CONTINUE
C
C ---------------------------------------------------------------------
C     PRINT ALL THE DIAGONALS OF  T  AND  R.
C ---------------------------------------------------------------------
  200 IF (MSG.LT.30) GO TO 260
      LENT = NROWRT*(NACTIV-1) + 1
      INCT = NROWRT - 1
      IF (NACTIV.GT.0) CALL DCOPY(NACTIV,RT(NACTIV,NCOLZ+1)
     *                             ,INCT,WRK1,1)
      IF (NACTIV.GT.0) THEN
         WRITE (REC,FMT=99988)
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         DO 220 I = 1, NACTIV, 5
            WRITE (REC,FMT=99994) (WRK1(J),J=I,MIN(I+4,NACTIV))
            CALL X04BAF(NOUT,REC(1))
  220    CONTINUE
      END IF
      IF (NCOLZ.GT.0) THEN
         WRITE (REC,FMT=99987)
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         DO 240 I = 1, NCOLZ, 5
            WRITE (REC,FMT=99994) (RT(J,J),J=I,MIN(I+4,NCOLZ))
            CALL X04BAF(NOUT,REC(1))
  240    CONTINUE
      END IF
C
  260 RETURN
C
C
C     END OF E04NAT  ( QPPRT )
99999 FORMAT (///' =================',/' QP ITERATION',I5,/' =========',
     *  '========')
99998 FORMAT (//'  ITN JDEL  JADD       STEP NHESS   OBJECTIVE NCOLZ N',
     *  'ORM GFREE  NORM ZTG   COND T COND ZHZ  HESS MOD')
99997 FORMAT (//'  ITN JDEL  JADD       STEP NHESS   OBJECTIVE NCOLZ  ',
     *  ' NORM QTG  NORM ZTG   COND T COND ZHZ  HESS MOD')
99996 FORMAT (I5,I5,A1,I5,A1,1P,D10.2,I6,1P,D12.4,I6,1P,D11.2,1P,D10.2,
     *  1P,2D9.1,1P,D10.2)
99995 FORMAT (/' QP VARIABLES')
99994 FORMAT (1P,5D15.6)
99993 FORMAT (/' STATUS OF THE QP BOUND CONSTRAINTS')
99992 FORMAT (1X,10I4)
99991 FORMAT (/' STATUS OF THE QP GENERAL CONSTRAINTS')
99990 FORMAT (/' LIST OF FREE QP VARIABLES')
99989 FORMAT (/' VALUES OF QP GENERAL LINEAR CONSTRAINTS')
99988 FORMAT (/' DIAGONALS OF QP WORKING SET FACTOR  T  ')
99987 FORMAT (/' DIAGONALS OF QP PRJ. HESSIAN FACTOR  R ')
      END
