      SUBROUTINE F01MAF(N,NZ,A,LA,INI,LINI,INJ,DROPTL,DENSW,W,IK,IW,
     *                  ABORT,INFORM,IFLAG)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     F01MAF PERFORMS AN INCOMPLETE CHOLESKY FACTORIZATION OF
C     A SPARSE SYMMETRIC POSITIVE DEFINITE MATRIX A.
C
C     THE ROUTINE IS BASED UPON THE HARWELL ROUTINE MA31A BY
C     N MUNKSGAARD.
C
C     THE DEVELOPMENT OF MA31A WAS SUPPORTED BY THE DANISH
C     NATURAL RESEARCH COUNCIL, GRANT NUMBER 511-10069
C
C     FOR A DESCRIPTION OF THE PARAMETERS AND USE OF THIS ROUTINE
C     SEE THE NAG LIBRARY MANUAL. PARAMETER ASSOCIATION IS AS BELOW.
C
C     ROUTINE    DOCUMENT
C     N           N
C     NZ         NZ
C     A           A
C     LA         LICN
C     INI        IRN
C     LINI       LIRN
C     INJ        ICN
C     DROPTL     DROPTL
C     DENSW      DENSW
C     INFORM     INFORM
C     W         WKEEP
C     IK        IKEEP
C     IW       IWORK
C     ABORT       ABORT
C     IFLAG       IFAIL
C
C     ARRAYS MA31I, MA31J AND MA31K CORRESPOND TO COMMON BLOCKS
C     OF THE SAME NAMES IN MA31A.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01MAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DENSW, DROPTL
      INTEGER           IFLAG, LA, LINI, N, NZ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LA), W(N,3)
      INTEGER           IK(N,2), INFORM(4), INI(LINI), INJ(LA), IW(N,6)
      LOGICAL           ABORT(3)
C     .. Local Scalars ..
      DOUBLE PRECISION  ONE, ZERO
      INTEGER           I, IAJ1, IERR, IFAIL, II, IR, J, K, KI, KJ, KK,
     *                  KL, KLL, KP, KPP, KR, LP, MP, NM1, NUAL, NZ0,
     *                  NZP1
      LOGICAL           L1, L2
C     .. Local Arrays ..
      DOUBLE PRECISION  MA31I(1)
      INTEGER           MA31J(5), MA31K(3)
      CHARACTER*1       P01REC(1)
      CHARACTER*80      REC(3)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01MAX, F01MAZ, X04AAF, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MOD, SQRT
C     .. Data statements ..
      DATA              ZERO/0.0D+0/, ONE/1.0D+0/
C     .. Executable Statements ..
C
      CALL X04AAF(0,LP)
      CALL X04ABF(0,MP)
C
      L1 = IFLAG/100 .EQ. 1
      L2 = MOD(IFLAG,100)/10 .EQ. 1
C
C     CHECK RESTRICTIONS ON INPUT PARAMETERS.
C
      IF (N.LT.1) GO TO 460
      IF (NZ.LT.N) GO TO 480
      IF (LINI.LT.NZ) GO TO 500
      IF (LA.LT.2*NZ) GO TO 520
      IF (ABS(DROPTL).GT.ONE) DROPTL = ONE
      MA31I(1) = DENSW
      IF (MA31I(1).GT.ONE .OR. MA31I(1).LT.ZERO) MA31I(1) = 0.8D0
C
C     INITIALIZE WORK ARRAYS.
C
      IERR = 0
      DO 40 I = 1, N
         DO 20 J = 1, 2
            W(I,J) = ZERO
            IK(I,J) = 0
   20    CONTINUE
   40 CONTINUE
      DO 60 I = 1, N
         W(I,3) = ZERO
         IW(I,5) = 0
   60 CONTINUE
      MA31K(3) = 0
      MA31J(4) = 0
      MA31J(3) = 0
C
C     COUNT NUMBER OF ELEMENTS
C
      DO 140 K = 1, NZ
         I = INI(K)
         J = INJ(K)
         IF (I.LT.1 .OR. I.GT.N) GO TO 540
         IF (J.LT.I .OR. J.GT.N) GO TO 540
         IF (I.EQ.J) GO TO 100
         IF (A(K).NE.ZERO) GO TO 80
         MA31J(4) = MA31J(4) + 1
         GO TO 140
   80    IK(I,1) = IK(I,1) + 1
         IK(J,2) = IK(J,2) + 1
         MA31K(3) = MA31K(3) + 1
         NUAL = MA31K(3)
         A(NUAL) = A(K)
         INI(NUAL) = I
         INJ(NUAL) = J
         GO TO 140
C
C        CHECK FOR DOUBLE ENTRIES ON THE DIAGONAL AND REMOVE
C        DIAGONAL FROM A TO W.
C
  100    MA31J(4) = MA31J(4) + 1
         IF (W(I,1).EQ.ZERO) GO TO 120
C
C        *****  ADVISORY MESSAGE
C
         IF (L1) THEN
            WRITE (REC,FMT=99999) I
            CALL X04BAF(MP,REC(1))
            CALL X04BAF(MP,REC(2))
            CALL X04BAF(MP,REC(3))
         END IF
         IF ( .NOT. ABORT(1)) GO TO 120
         IERR = P01ABF(IFLAG,5,SRNAME,0,P01REC)
         GO TO 620
  120    W(I,1) = W(I,1) + A(K)
         IF (W(I,1).LE.ZERO) GO TO 560
  140 CONTINUE
C
C     NZ0 IS THE NUMBER OF OFF DIAGONAL NON - ZEROS.
C
      NZ0 = NZ - MA31J(4)
      MA31J(2) = NZ0
      MA31J(1) = NZ0
C
C     ON DIAGONAL MATRIX MAKE SPECIAL EXIT.
C
      IF (NZ0.EQ.0) GO TO 660
C
C     INITIALIZE IW( I, 1 )AND IW( I, 2 )TO POINT JUST BEYOND WHERE
C     THE LAST COMPONENT OF ROW/COLUMN I WILL BE STORED.
C
      KJ = LINI - NZ0 + 1
      KI = 1
      DO 160 I = 1, N
         KI = KI + IK(I,1)
         KJ = KJ + IK(I,2)
         IW(I,1) = KI
         IW(I,2) = KJ
         W(I,1) = ONE/SQRT(W(I,1))
  160 CONTINUE
C
C     REORDER BY ROWS USING IN - PLACE SORT ALGORITHM.
C
      CALL F01MAX(INI,INJ,NZ0,IW,N,A)
C
C     CHECK FOR DOUBLE ENTRIES WHILE USING THE CONSTRUCTED ROW FILE
C     TO SET UP THE COLUMN FILE AND COMPRESS THE ROWFILE.
C
      KK = 0
      DO 260 IR = 1, N
         KPP = IW(IR,1)
         IW(IR,1) = KK + 1
         KLL = KPP + IK(IR,1) - 1
         IF (KPP.GT.KLL) GO TO 260
C
C        LOAD ROW IR INTO W( *, 3 ).
C
         DO 200 K = KPP, KLL
            J = INJ(K)
            IF (W(J,3).EQ.ZERO) GO TO 180
C
C           *****  ADVISORY MESSAGE
C
            IF (L1) THEN
               WRITE (REC,FMT=99998) IR, J
               CALL X04BAF(MP,REC(1))
               CALL X04BAF(MP,REC(2))
               CALL X04BAF(MP,REC(3))
            END IF
            IF ( .NOT. ABORT(1)) GO TO 180
            IERR = P01ABF(IFLAG,5,SRNAME,0,P01REC)
            GO TO 620
  180       W(J,3) = W(J,3) + A(K)
  200    CONTINUE
C
C        RELOAD ROW IR INTO ARRAYS A AND INJ AND ADJUST INI.
C
         DO 240 K = KPP, KLL
            J = INJ(K)
            IF (W(J,3).EQ.ZERO) GO TO 220
            KK = KK + 1
            A(KK) = W(J,3)*W(IR,1)*W(J,1)
            INJ(KK) = J
            W(J,3) = ZERO
            KR = IW(J,2) - 1
            IW(J,2) = KR
            INI(KR) = IR
            GO TO 240
  220       MA31J(4) = MA31J(4) + 1
            MA31J(1) = MA31J(1) - 1
            MA31J(2) = MA31J(2) - 1
            IK(IR,1) = IK(IR,1) - 1
            IK(J,2) = IK(J,2) - 1
  240    CONTINUE
  260 CONTINUE
      IF (IERR.NE.5) GO TO 320
      NZ0 = NZ - MA31J(4)
C
C     ZERO UNUSED LOCATIONS IN INI.
C
      NM1 = N - 1
      DO 300 I = 1, NM1
         KP = IW(I,2) + IK(I,2)
         KL = IW(I+1,2) - 1
         IF (KP.GT.KL) GO TO 300
         DO 280 K = KP, KL
            INI(K) = 0
  280    CONTINUE
  300 CONTINUE
C
C     STORE INPUT MATRIX.
C
  320 NUAL = LA + 1
      DO 380 II = 1, N
         IW(II,6) = IK(II,1)
         I = N - II + 1
         W(I,2) = ONE
         KP = IW(I,1)
         KL = KP + IK(I,1) - 1
         IF (KP.GT.KL) GO TO 360
         DO 340 KK = KP, KL
            K = KP + KL - KK
            NUAL = NUAL - 1
            A(NUAL) = A(K)
            INJ(NUAL) = INJ(K)
  340    CONTINUE
  360    IW(I,1) = NUAL - NZ0
  380 CONTINUE
C
C     SET DIFFERENT PARAMETERS.
C
      MA31K(1) = 0
      NZP1 = NZ0 + 1
      MA31K(2) = IW(1,2)
      MA31K(3) = NUAL - NZ0
C
C     ACTIVATE INCOMPLETE FACTORIZATION.
C
      IAJ1 = LA - NZ0
      IFAIL = IFLAG
      CALL F01MAZ(N,NZ0,W(1,2),A(NZP1),INI,INJ(NZP1)
     *            ,LINI,IAJ1,IK,IW,IW(1,3),W(1,3)
     *            ,DROPTL,MA31I,MA31J,MA31K,ABORT,IFLAG)
      IF (MA31J(1).GT.IAJ1 .OR. MA31J(2).GT.LINI) GO TO 580
      IF (IFLAG.EQ.6 .OR. IFLAG.EQ.7) GO TO 640
      IERR = IFLAG + IERR
      IFLAG = IFAIL
C
C     THE FACTORIZATION IS TERMINATED.
C
      KP = 1
      DO 440 I = 1, N
         KL = KP + IW(I,6) - 1
         IF (KP.GT.KL) GO TO 420
         DO 400 K = KP, KL
            INI(K) = I
  400    CONTINUE
  420    KP = KL + 1
  440 CONTINUE
      GO TO 620
C
C     THE FOLLOWING INSTRUCTIONS IMPLEMENT THE FAILURE EXITS.
C
C     *****  ERROR MESSAGES
C
  460 IF (L2) THEN
         WRITE (REC,FMT=99991)
         CALL X04BAF(LP,REC(1))
         WRITE (REC,FMT=99997)
         CALL X04BAF(LP,REC(1))
         CALL X04BAF(LP,REC(2))
         CALL X04BAF(LP,REC(3))
      END IF
      IERR = 1
      GO TO 600
  480 IF (L2) THEN
         WRITE (REC,FMT=99991)
         CALL X04BAF(LP,REC(1))
         WRITE (REC,FMT=99996)
         CALL X04BAF(LP,REC(1))
         CALL X04BAF(LP,REC(2))
         CALL X04BAF(LP,REC(3))
      END IF
      IERR = 1
      GO TO 600
  500 IF (L2) THEN
         WRITE (REC,FMT=99991)
         CALL X04BAF(LP,REC(1))
         WRITE (REC,FMT=99995)
         CALL X04BAF(LP,REC(1))
         CALL X04BAF(LP,REC(2))
         CALL X04BAF(LP,REC(3))
      END IF
      IERR = 1
      GO TO 600
  520 IF (L2) THEN
         WRITE (REC,FMT=99991)
         CALL X04BAF(LP,REC(1))
         WRITE (REC,FMT=99994)
         CALL X04BAF(LP,REC(1))
         CALL X04BAF(LP,REC(2))
         CALL X04BAF(LP,REC(3))
      END IF
      IERR = 1
      GO TO 600
  540 IF (L2) THEN
         WRITE (REC,FMT=99991)
         CALL X04BAF(LP,REC(1))
         WRITE (REC,FMT=99993) K, I, J
         CALL X04BAF(LP,REC(1))
         CALL X04BAF(LP,REC(2))
         CALL X04BAF(LP,REC(3))
      END IF
      IERR = 2
      GO TO 600
  560 IF (L2) THEN
         WRITE (REC,FMT=99991)
         CALL X04BAF(LP,REC(1))
         WRITE (REC,FMT=99992) I
         CALL X04BAF(LP,REC(1))
         CALL X04BAF(LP,REC(2))
         CALL X04BAF(LP,REC(3))
      END IF
      IERR = 3
      GO TO 600
  580 IF (L2) THEN
         WRITE (REC,FMT=99991)
         CALL X04BAF(LP,REC(1))
         WRITE (REC,FMT=99990)
         CALL X04BAF(LP,REC(1))
         CALL X04BAF(LP,REC(2))
         CALL X04BAF(LP,REC(3))
      END IF
      IERR = 4
      IFLAG = IFAIL
  600 IERR = P01ABF(IFLAG,IERR,SRNAME,0,P01REC)
C
  620 IFLAG = IERR
  640 INFORM(1) = MA31J(1)
      INFORM(2) = MA31J(2)
      INFORM(3) = MA31J(4)
      INFORM(4) = MA31J(5)
      RETURN
C
C     THE REMAINING PART OF THE CODE IMPLEMENTS THE DIAGONAL MATRIX
C     CASE.
C
  660 DO 680 I = 1, N
         IK(I,2) = I
         IK(I,1) = 0
         IW(I,5) = 0
         W(I,1) = ONE/SQRT(W(I,1))
         W(I,2) = ONE
  680 CONTINUE
      GO TO 620
C
C
C     END OF F01MAF.
C
99999 FORMAT (//' WARNING MORE THAN ONE DIAGONAL ENTRY IN ROW',I5)
99998 FORMAT (//' WARNING THERE IS MORE THAN ONE ENTRY IN ROW',I5,' AN',
     *  'D COLUMN',I5)
99997 FORMAT (//34X,' N  IS OUT OF RANGE.')
99996 FORMAT (//34X,' NZ IS OUT OF RANGE.')
99995 FORMAT (//34X,'LINDI IS OUT OF RANGE.')
99994 FORMAT (//34X,'LA IS OUT OF RANGE.')
99993 FORMAT (//34X,'ELEMENT',I7,' IS IN ROW',I5,' AND COLUMN',I5)
99992 FORMAT (//34X,'DIAGONAL ELEMENT',I5,' IS ZERO OR NEGATIVE')
99991 FORMAT (' ERROR RETURN FROM F01MAF BECAUSE')
99990 FORMAT (//34X,'NO FURTHER FILL IN CAN BE ALLOWED')
      END
