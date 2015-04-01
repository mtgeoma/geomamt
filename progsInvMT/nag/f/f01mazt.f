      SUBROUTINE F01MAZ(N,NZ,D,A,INI,INJ,IAI,IAJ,IK,IP,IW,W,C,MA31I,
     *                  MA31J,MA31K,ABORT,IFLAG)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 17 REVISED. IER-1706. (OCT 1995).
C
C     IP( I, 1 ), IP( I, 2 )POINT TO THE START OF ROW/COLUMN I.
C     IK( I, 1 ), IK( I, 2 )HOLD THE NUMBER OF NONZEROS IN ROW/COLUMN
C     I OF THE LOWER TRIANGULAR PART OF A.
C     DURING THE BODY OF THIS ROUTINE THE VECTORS IW( *, 3 ), IW( *, 1 )
C     IW( *, 2 )ARE USED TO HOLD DOUBLY LINKED LISTS OF ROWS THAT HAVE
C     NOT BEEN PIVOTAL AND HAVE EQUAL NUMBER OF NONZEROS.
C     IW( I, 3 )HOLD FIRST ROW/COLUMN TO HAVE I NONZEROS OR ZERO IF
C     THERE ARE NONE.
C     IW( I, 1 )HOLD ROW/COLUMN NUMBER OF ROW/COLUMN PRIOR TO ROW I IN
C     ITS LIST OR ZERO IF NONE.
C     IW( I, 2 )HOLD ROW/COLUMN NUMBER OF ROW/COLUMN AFTER ROW I IN ITS
C     LIST OR ZERO IF NONE.
C     DURING THE MAIN BODY OF THE SUBROUTINE INI AND INJ KEEP A COLUMN
C     FILE AND A ROW FILE CONTAINING RESPECTIVELY THE ROW NUMBERS OF
C     THE NON - ZEROS OF EACH COLUMN AND THE COLUMN NUMBERS OF THE  NON
C     ZEROS OF EACH ROW. THE IP ARRAYS POINT TO THE START POSITION IN
C     INI AND INJ OF EACH COLUMN AND ROW.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01MAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  C
      INTEGER           IAI, IAJ, IFLAG, N, NZ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IAJ), D(N), MA31I(1), W(N)
      INTEGER           IK(N,2), INI(IAI), INJ(IAJ), IP(N,2), IW(N,3),
     *                  MA31J(5), MA31K(3)
      LOGICAL           ABORT(3)
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION  AA, AL, ALFA, B1, B2, C0, CMAX, DJP, EPSTOL,
     *                  ONE, PFILL, PIVT, ZERO
      INTEGER           I, IERR, IFAIL, II, IIP, IIPP1, IL, ILOCS, IN,
     *                  IP1, IPDP1, IR, J, J1, JJ, JP, K, KC, KK, KL,
     *                  KL1, KL2, KLC, KLL, KLR, KP, KP2, KPC, KPI, KPJ,
     *                  KPP, KPR, KR, KRL, KS, L, LFULDD, LFULL, MCL,
     *                  MP, NC, NFILL, NLOCS, NM1, NR, NRJP, NRM1, NUAL,
     *                  NUCL, NURL, NZ0, NZC, NZI
      LOGICAL           CHANGE, L1
C     .. Local Arrays ..
      CHARACTER         P01REC(1)
      CHARACTER*80      REC(4)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, F01MAY, X02ZAZ, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, INT, MAX, SQRT
C     .. Common blocks ..
      COMMON            /AX02ZA/WMACH
C     .. Save statements ..
      SAVE              /AX02ZA/
C     .. Data statements ..
      DATA              ZERO/0.0D+0/, ONE/1.0D+0/
C
C     .. Executable Statements ..
      CALL X02ZAZ
C
      CMAX = WMACH(5)
      EPSTOL = 100.0D0*WMACH(3)
C
C     INITIALIZE IW( *, 3 )AND LOCAL VARIABLES.
C
      CALL X04ABF(0,MP)
C
      IFAIL = IFLAG
      IFLAG = 0
      L1 = IFAIL .GE. 100
      CHANGE = .TRUE.
      IF (C.LE.ZERO) CHANGE = .FALSE.
      NZ0 = NZ
      MA31J(5) = N
      ALFA = 1.0D0/0.90D0
      B1 = -0.03D0
      B2 = 0.03D0
      NFILL = IAJ - NZ0 - N
      MCL = MA31J(2)
      C0 = ABS(C)
      C0 = C0/DBLE(N)
      C = C**2
      DO 20 I = 1, N
         D(I) = D(I) + C0
         IW(I,3) = 0
   20 CONTINUE
      NLOCS = NZ
C
C     SET UP  LINKED LISTS OF ROWS/COLUMNS WITH EQUAL NUMBER OF
C     NON - ZEROS.
C
      DO 40 I = 1, N
         NZI = IK(I,1) + IK(I,2) + 1
         IN = IW(NZI,3)
         IW(NZI,3) = I
         IW(I,2) = IN
         IW(I,1) = 0
         IF (IN.NE.0) IW(IN,1) = I
   40 CONTINUE
C
C     START THE ELIMINATION LOOP
C
      DO 1180 IIP = 1, N
C
C        SEARCH ROWS WITH NRJP NONZEROS.
C
         DO 60 NRJP = 1, N
            JP = IW(NRJP,3)
            IF (JP.GT.0) GO TO 80
   60    CONTINUE
C
C        ROW JP  IS USED AS PIVOT.
C
C        REMOVE ROWS/COLUMNS INVOLVED IN ELIMINATION FROM ORDERING
C        VECTORS.
C
   80    CONTINUE
         DO 200 L = 1, 2
            KPP = IP(JP,L)
            KLL = IK(JP,L) + KPP - 1
            IF (KPP.GT.KLL) GO TO 200
            DO 180 K = KPP, KLL
               IF (L.EQ.2) GO TO 100
               J = INJ(K)
               GO TO 120
  100          CONTINUE
               J = INI(K)
  120          CONTINUE
               IL = IW(J,1)
               IN = IW(J,2)
               IW(J,2) = -1
               IF (IN.LT.0) GO TO 180
               IF (IL.EQ.0) GO TO 140
               IW(IL,2) = IN
               GO TO 160
  140          CONTINUE
               NZ = IK(J,1) + IK(J,2) + 1
               IW(NZ,3) = IN
  160          CONTINUE
               IF (IN.GT.0) IW(IN,1) = IL
  180       CONTINUE
  200    CONTINUE
C
C        REMOVE JP     FROM ORDERING VECTORS
C
         IL = IW(JP,1)
         IN = IW(JP,2)
         IW(JP,2) = -10
         IF (IN.LT.0) GO TO 220
         NZ = IK(JP,1) + IK(JP,2) + 1
         IW(NZ,3) = IN
         IF (IN.GT.0) IW(IN,1) = IL
C
C        STORE  PIVOT.
C
  220    CONTINUE
         IW(JP,1) = -IIP
C
C        COMPRESS ROW FILE IF NECESSARY.
C
         IF (MA31J(1)+IK(JP,1)+IK(JP,2).GT.IAJ-N) C = CMAX
         IF ((MA31K(1)+IK(JP,1)+IK(JP,2)).LT.MA31K(3)) GO TO 260
         CALL F01MAY(A,INJ,IAJ,N,IK,IP,.TRUE.,MA31J,MA31K)
         NLOCS = 0
         DO 240 ILOCS = 1, N
            NLOCS = NLOCS + IK(ILOCS,1)
  240    CONTINUE
  260    CONTINUE
         KP = IP(JP,1)
         KL = IK(JP,1) + KP - 1
         IP(JP,1) = MA31K(1) + 1
         NURL = MA31K(1)
         IF (KP.GT.KL) GO TO 380
C
C        REMOVE JP FROM COLUMNS CONTAINED IN THE PIVOT ROW.
C
         DO 360 K = KP, KL
            J = INJ(K)
            KPC = IP(J,2)
            NZ = IK(J,2) - 1
            IK(J,2) = NZ
            KLC = KPC + NZ
            IF (KLC.GT.KPC) GO TO 280
            INI(KPC) = 0
            GO TO 340
  280       CONTINUE
            DO 300 KC = KPC, KLC
               IF (JP.EQ.INI(KC)) GO TO 320
  300       CONTINUE
            KC = KLC + 1
  320       CONTINUE
            INI(KC) = INI(KLC)
            INI(KLC) = 0
  340       CONTINUE
            MA31J(2) = MA31J(2) - 1
            NURL = NURL + 1
            INJ(NURL) = J
            A(NURL) = A(K)
            INJ(K) = 0
  360    CONTINUE
C
C        TRANSFORM COLUMN PART OF PIVOT ROW TO THE ROW FILE.
C
  380    CONTINUE
         KP2 = IP(JP,2)
         KL2 = IK(JP,2) + KP2 - 1
         IF (KP2.GT.KL2) GO TO 460
         DO 440 K = KP2, KL2
            NURL = NURL + 1
            MA31J(2) = MA31J(2) - 1
            I = INI(K)
            KPR = IP(I,1)
            KLR = KPR + IK(I,1) - 1
            DO 400 KR = KPR, KLR
               IF (JP.EQ.INJ(KR)) GO TO 420
  400       CONTINUE
            KR = KLR + 1
  420       CONTINUE
            INJ(KR) = INJ(KLR)
            A(NURL) = A(KR)
            A(KR) = A(KLR)
            INJ(KLR) = 0
            IK(I,1) = IK(I,1) - 1
            NLOCS = NLOCS - 1
            INJ(NURL) = I
            INI(K) = 0
  440    CONTINUE
  460    CONTINUE
         NZC = IK(JP,1) + IK(JP,2)
         MA31K(1) = NURL
         NLOCS = NLOCS + NZC - IK(JP,1)
         IK(JP,1) = NZC
         IK(JP,2) = 0
C
C        UNPACK PIVOT ROW AND CONTROL DIAGONAL VALUE.
C
         KP = IP(JP,1)
         KL = KP + NZC - 1
         C0 = EPSTOL
         IF (KP.GT.KL) GO TO 500
         DO 480 K = KP, KL
            AA = A(K)
            C0 = MAX(C0,ABS(AA))
            J = INJ(K)
            W(J) = AA
  480    CONTINUE
  500    CONTINUE
         DJP = D(JP)
         IF (DJP.GT.C0/100.0D0) GO TO 540
C
C        *****  ADVISORY MESSAGE
C
         IF (L1) THEN
            WRITE (REC,FMT=99999) JP
            CALL X04BAF(MP,REC(1))
            CALL X04BAF(MP,REC(2))
            CALL X04BAF(MP,REC(3))
            CALL X04BAF(MP,REC(4))
         END IF
         IF ( .NOT. ABORT(2)) GO TO 520
         IERR = P01ABF(IFAIL,6,SRNAME,0,P01REC)
         IFLAG = IERR
         RETURN
  520    CONTINUE
         D(JP) = C0
         IF (C0.EQ.EPSTOL) D(JP) = ONE
  540    CONTINUE
         IF (KP.GT.KL) GO TO 1160
C
C        PERFORM ROW OPERATIONS.
C
         DO 1080 NC = 1, NZC
            KC = IP(JP,1) + NC - 1
            IR = INJ(KC)
            AL = A(KC)/D(JP)
C
C           COMPRESS ROW FILE IF NECESSARY.
C
            IF (MA31J(1)+IK(IR,1)+IK(JP,1).GT.IAJ-N) C = CMAX
            IF ((MA31K(1)+IK(IR,1)+IK(JP,1)).LT.MA31K(3)) GO TO 580
            CALL F01MAY(A,INJ,IAJ,N,IK,IP,.TRUE.,MA31J,MA31K)
            NLOCS = 0
            DO 560 ILOCS = 1, N
               NLOCS = NLOCS + IK(ILOCS,1)
  560       CONTINUE
  580       CONTINUE
            KR = IP(IR,1)
            KRL = KR + IK(IR,1) - 1
            IF (KR.GT.KRL) GO TO 620
C
C           SCAN THE OTHER ROW AND CHANGE SIGN IN IW FOR EACH
C           COMMON COLUMN NUMBER
C
            DO 600 KS = KR, KRL
               J = INJ(KS)
               IF (IW(J,2).NE.(-1)) GO TO 600
               IW(J,2) = 1
               A(KS) = A(KS) - AL*W(J)
  600       CONTINUE
C
C           SCAN PIVOT ROW FOR FILLS.
C
  620       CONTINUE
            DO 1060 KS = KP, KL
               J = INJ(KS)
C
C              ONLY ENTRIES IN THE UPPER TRIANGULAR PART ARE CONSIDERED.
C
               IF (J.LT.IR) GO TO 1040
               IF (IW(J,2).EQ.1) GO TO 1040
               AA = -AL*W(J)
               IF (IR.NE.J) GO TO 640
               D(IR) = D(IR) + AA
               GO TO 1040
  640          CONTINUE
               IF (AA*AA.GT.C*ABS(D(IR)*D(J))) GO TO 660
               D(J) = D(J) + AA
               D(IR) = D(IR) + AA
               GO TO 1040
  660          CONTINUE
               MA31J(1) = MA31J(1) + 1
               IF (MA31J(1).LE.IAJ) GO TO 680
               RETURN
  680          CONTINUE
               IK(IR,1) = IK(IR,1) + 1
               NLOCS = NLOCS + 1
C
C              IF POSSIBLE PLACE THE NEW ELEMENT NEXT TO THE PRESENT
C              ENTRY.
C
C              TRY IF THERE IS ROOM AT THE END OF THE ENTRY.
C
               IF (KR.GT.KRL) GO TO 800
               IF (KRL.EQ.IAJ) GO TO 700
               IF (INJ(KRL+1).NE.0) GO TO 700
               KRL = KRL + 1
               INJ(KRL) = J
               A(KRL) = AA
               GO TO 820
C
C              TRY IF THERE IS ROOM AHEAD OF PRESENT ENTRY.
C
  700          CONTINUE
               IF (KR.NE.MA31K(3)) GO TO 720
               MA31K(3) = MA31K(3) - 1
               GO TO 740
  720          CONTINUE
               IF (INJ(KR-1).NE.0) GO TO 760
  740          CONTINUE
               KR = KR - 1
               IP(IR,1) = KR
               INJ(KR) = J
               A(KR) = AA
               GO TO 820
C
C              NEW ENTRY HAS TO BE CREATED.
C
  760          CONTINUE
               NUAL = MA31K(3)
               DO 780 KK = KR, KRL
                  NUAL = NUAL - 1
                  INJ(NUAL) = INJ(KK)
                  A(NUAL) = A(KK)
                  INJ(KK) = 0
  780          CONTINUE
               MA31K(3) = NUAL
C
C              ADD THE NEW ELEMENT.
C
  800          CONTINUE
               MA31K(3) = MA31K(3) - 1
               NUAL = MA31K(3)
               INJ(NUAL) = J
               A(NUAL) = AA
               IP(IR,1) = NUAL
               KR = NUAL
               KRL = KR + IK(IR,1) - 1
C
C              IF ARRAYS IP(I,1) AND IK(I,1) ARE FULL ,
C              IE THERE IS NO ROOM FOR THE
C              NEW ENTRY , SET MA31J(1) GREATER THAN IAJ AND RETURN.
C              IP(I,1) POINTS TO STARTING POSITION OF ELEMENTS IN ROW I.
C              IK(I,1) CONTAINS THE NUMBER OF ELEMENTS IN ROW I.
C              NOTE - LABEL 740 WAS MOVED FROM LINE CONTAINING
C              NZ = IK(J,2)
C
  820          CONTINUE
               IF (NLOCS.GT.IAJ) THEN
                  MA31J(1) = IAJ + 1
                  RETURN
               END IF
C
C              CREATE FILL IN COLUMN FILE.
C
               NZ = IK(J,2)
               K = IP(J,2)
               KL1 = K + NZ - 1
               MA31J(2) = MA31J(2) + 1
               IF (MA31J(2).LE.IAI) GO TO 840
               RETURN
C
C              IF POSSIBLE PLACE NEW ELEMENT AT THE END OF
C              PRESENT ENTRY.
C
  840          CONTINUE
               IF (NZ.EQ.0) GO TO 920
               IF (KL1.EQ.IAI) GO TO 860
               IF (INI(KL1+1).NE.0) GO TO 860
               INI(KL1+1) = IR
               GO TO 1020
C
C              IF POSSIBLE PLACE ELEMENT AHEAD OF PRESENT ENTRY.
C
  860          CONTINUE
               IF (K.NE.MA31K(2)) GO TO 880
               IF (MA31K(2).EQ.1) GO TO 920
               MA31K(2) = MA31K(2) - 1
               GO TO 900
  880          CONTINUE
               IF (INI(K-1).NE.0) GO TO 920
  900          CONTINUE
               K = K - 1
               INI(K) = IR
               IP(J,2) = K
               GO TO 1020
C
C              NEW ENTRY HAS TO BE CREATED.
C
  920          CONTINUE
               IF (NZ+1.LT.MA31K(2)) GO TO 940
C
C              COMPRESS COLUMN FILE IF THERE IS NOT ROOM FOR NEW ENTRY.
C
               IF (MA31J(2)+NZ+2.GE.IAI) C = CMAX
               CALL F01MAY(A,INI,IAI,N,IK(1,2),IP(1,2),.FALSE.,MA31J,
     *                     MA31K)
               K = IP(J,2)
               KL1 = K + NZ - 1
C
C              TRANSFER OLD ENTRY INTO NEW.
C
  940          CONTINUE
               IF (K.GT.KL1) GO TO 1000
               NUCL = MA31K(2)
               IF (KL1-K.NE.NUCL-1) GO TO 960
               MA31J(2) = IAI + 1
               RETURN
  960          CONTINUE
               DO 980 KK = K, KL1
                  NUCL = NUCL - 1
                  INI(NUCL) = INI(KK)
                  INI(KK) = 0
  980          CONTINUE
               MA31K(2) = NUCL
C
C              ADD THE NEW ELEMENT.
C
 1000          CONTINUE
               MA31K(2) = MA31K(2) - 1
               NUCL = MA31K(2)
               INI(NUCL) = IR
               IP(J,2) = NUCL
 1020          CONTINUE
               IK(J,2) = NZ + 1
 1040          CONTINUE
               IW(J,2) = -1
 1060       CONTINUE
 1080    CONTINUE
C
C        UPDATE ORDERING ARRAYS.
C
         DO 1100 K = KP, KL
            J = INJ(K)
            W(J) = 0.D0
            A(K) = A(K)/D(JP)
            NZ = IK(J,1) + IK(J,2) + 1
            IN = IW(NZ,3)
            IW(J,2) = IN
            IW(J,1) = 0
            IW(NZ,3) = J
            IF (IN.NE.0) IW(IN,1) = J
 1100    CONTINUE
         MCL = MAX(MCL,MA31J(2))
         PIVT = DBLE(IIP)/DBLE(N)
C
C        GIVE WARNING IF AVAILABLE SPACE IS USED TOO EARLY.
C
         IF (C.NE.CMAX) GO TO 1120
         IF (MA31J(5).LT.IIP) GO TO 1160
         MA31J(5) = IIP
         IF (PIVT.GT.0.9D0) GO TO 1160
C
C        ***** ADVISORY MESSAGE
C
         IF (L1) THEN
            WRITE (REC,FMT=99998) IIP
            CALL X04BAF(MP,REC(1))
            CALL X04BAF(MP,REC(2))
            CALL X04BAF(MP,REC(3))
         END IF
         IF ( .NOT. ABORT(3)) GO TO 1160
         IERR = P01ABF(IFAIL,7,SRNAME,0,P01REC)
         IFLAG = IERR
         RETURN
C
C        CHANGE C IF NECESSARY.
C
 1120    CONTINUE
         IF ( .NOT. CHANGE) GO TO 1160
         PFILL = DBLE(MA31J(1)-NZ0)/DBLE(NFILL)
         IF (PIVT.GE.0.9D0) GO TO 1160
         IF (PFILL.LT.ALFA*PIVT+B1) GO TO 1140
         IF (PFILL.LT.ALFA*PIVT+B2) GO TO 1160
         C = 2.25D0*C
 1140    CONTINUE
         ALFA = (1.0D0-PFILL)/(0.9D0-PIVT)
         B1 = PFILL - PIVT*ALFA - 0.03D0
         B2 = B1 + 0.06D0
C
C        IF THE MATRIX IS FULL  THEN STOP THE SPARSE ANALYZE.
C
 1160    CONTINUE
         NR = N - IIP
         LFULL = (NR*(NR-1)/2)
         LFULDD = INT(MA31I(1)*DBLE(LFULL))
         IF (MA31J(2).GE.LFULDD .AND. MA31K(1)+LFULL.LT.IAJ) GO TO 1200
 1180 CONTINUE
C
C     ELIMINATION LOOP TERMINATES
C     AFTER DEVIATION WE FACTORIZE THE REMAINING FULL MATRIX.
C
 1200 CONTINUE
      MA31J(5) = IIP
      C = SQRT(C)
      MA31J(2) = MCL
      IF ( .NOT. CHANGE) C = -C
C
C     THE ORDER OF THE FULL MATRIX IS NR.
C     LOOP THROUGH ROWS IN THE ACTIVE MATRIX AND STORE ROW
C     NUMBERS IN INI.
C
      KK = 0
      DO 1260 I = 1, NR
         JP = IW(I,3)
 1220    CONTINUE
         IF (JP.LE.0) GO TO 1240
         KK = KK + 1
         INI(KK) = JP
         JP = IW(JP,2)
         GO TO 1220
 1240    CONTINUE
         IF (KK.EQ.NR) GO TO 1280
 1260 CONTINUE
C
C     MAKE A SORT OF ROW NUMBERS IN INI.
C
 1280 CONTINUE
      IF (NR.EQ.1) GO TO 1340
      NRM1 = NR - 1
      DO 1320 I = 1, NRM1
         J1 = I + 1
         DO 1300 J = J1, NR
            IF (INI(J).GT.INI(I)) GO TO 1300
            JJ = INI(I)
            INI(I) = INI(J)
            INI(J) = JJ
 1300    CONTINUE
 1320 CONTINUE
 1340 CONTINUE
      DO 1360 I = 1, NR
         II = INI(I)
         IW(II,1) = -(MA31J(5)+I)
 1360 CONTINUE
C
C     MAKE AN ORDERED LIST OF THE PIVOTS.
C
      DO 1380 I = 1, N
         IR = -IW(I,1)
         IK(IR,2) = I
 1380 CONTINUE
C
C     MOVE FULL MATRIX TO THE FRONT AND ORDER.
C
      IPDP1 = MA31J(5) + 1
      NM1 = N - 1
      IF (IPDP1.GT.NM1) RETURN
      DO 1540 IIP = IPDP1, NM1
         JP = IK(IIP,2)
         KP = IP(JP,1)
         KL = KP + IK(JP,1) - 1
C
C        MOVE ROW JP TO W.
C
         IF (KP.GT.KL) GO TO 1420
         DO 1400 K = KP, KL
            AA = A(K)
            C0 = MAX(C0,ABS(AA))
            J = INJ(K)
            INJ(K) = 0
            W(J) = A(K)
 1400    CONTINUE
 1420    CONTINUE
         DJP = D(JP)
         IF (DJP.GT.C0/100.0D0) GO TO 1460
C        ***** ADVISORY MESSAGE.
         IF (L1) THEN
            WRITE (REC,FMT=99999) JP
            CALL X04BAF(MP,REC(1))
            CALL X04BAF(MP,REC(2))
            CALL X04BAF(MP,REC(3))
            CALL X04BAF(MP,REC(4))
         END IF
         IF ( .NOT. ABORT(2)) GO TO 1440
         IERR = P01ABF(IFAIL,6,SRNAME,0,P01REC)
         IFLAG = IERR
         RETURN
 1440    CONTINUE
         D(JP) = C0
         IF (C0.EQ.EPSTOL) D(JP) = ONE
C
C        COMPRESS FILE IF NECESSARY.
C
 1460    CONTINUE
         IF (MA31K(1)+N-IIP.LT.MA31K(3)) GO TO 1500
         CALL F01MAY(A,INJ,IAJ,N,IK,IP,.TRUE.,MA31J,MA31K)
         NLOCS = 0
         DO 1480 ILOCS = 1, N
            NLOCS = NLOCS + IK(ILOCS,1)
 1480    CONTINUE
 1500    CONTINUE
         IP(JP,1) = MA31K(1) + 1
         NLOCS = NLOCS + N - IIP - IK(JP,1)
         IK(JP,1) = N - IIP
C
C        MOVE ROWS AND COLUMN INDICES INTO PIVOTAL ORDER.
C
         NURL = MA31K(1)
         IIPP1 = IIP + 1
         DO 1520 I = IIPP1, N
            J = IK(I,2)
            NURL = NURL + 1
            A(NURL) = W(J)
            INJ(NURL) = J
            W(J) = ZERO
 1520    CONTINUE
         MA31K(1) = NURL
 1540 CONTINUE
      MA31J(1) = MA31K(1)
C
C     FACTORIZE THE FULL MATRIX.
C
      DO 1600 IIP = IPDP1, NM1
         JP = IK(IIP,2)
         KPI = IP(JP,1)
         IP1 = IIP + 1
         IF (IP1.EQ.N) GO TO 1580
C
C        LOOP THROUGH THE OTHER ROW
C
         DO 1560 J = IP1, NM1
            JJ = IK(J,2)
            KPJ = IP(JJ,1)
            AL = A(KPI)/D(JP)
C
C           BEWARE THE FOLLOWING STATEMENT IF PASS BY COPY EVER COMES.
C
            CALL DAXPY(IK(JJ,1),-AL,A(KPI+1),1,A(KPJ),1)
C
            D(JJ) = D(JJ) - AL*A(KPI)
C
C           STORE FACTOR AND PROCEED TO NEXT ROW.
C
            A(KPI) = AL
            KPI = KPI + 1
 1560    CONTINUE
C
C        MODIFY LAST DIAGONAL ENTRY
C
 1580    CONTINUE
         JJ = IK(N,2)
         AL = A(KPI)/D(JP)
         D(JJ) = D(JJ) - AL*A(KPI)
         IF (D(JJ).LE.C0/100.0D0) GO TO 1620
         A(KPI) = AL
 1600 CONTINUE
      RETURN
C     ***** ADVISORY MESSAGE.
 1620 CONTINUE
      IF (L1) THEN
         WRITE (REC,FMT=99999) JP
         CALL X04BAF(MP,REC(1))
         CALL X04BAF(MP,REC(2))
         CALL X04BAF(MP,REC(3))
         CALL X04BAF(MP,REC(4))
      END IF
      IF ( .NOT. ABORT(2)) GO TO 1640
      IERR = P01ABF(IFAIL,6,SRNAME,0,P01REC)
      IFLAG = IERR
      RETURN
 1640 CONTINUE
      D(JJ) = C0
      IF (C0.EQ.EPSTOL) D(JJ) = ONE
      A(KPI) = AL
      RETURN
C
C
C     END OF F01MAZ. ( MA31C. )
C
99999 FORMAT (//'  WARNING   MODIFICATION OF ZERO OR NEGATIVE DIAGONAL',
     *       ' ENTRY HAS BEEN PERFORMED',/12X,'IN LOCATION',I7)
99998 FORMAT (//' WARNING  AVAILABLE SPACE USED AT PIVOT STEP',I7)
      END
