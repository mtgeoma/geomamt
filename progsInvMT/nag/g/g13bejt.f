      SUBROUTINE G13BEJ(XXY,IXXY,N,NXSP,PXS,IPXS,WDS,IWDS,NFR,MT,W,IDW,
     *                  STTF,ISTTF,NSTTF,MR,C,KZEF,BF,NBFQ,NBF,KSS,IERR)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     SUBROUTINE G13BEJ CONSTRUCTS THE T.F. TYPE STATE SET AND
C     ALLOWS OPTIONAL REPLACEMENT OF X(T) AND Y(T) BY Z(T) AND E(T)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  C
      INTEGER           IDW, IERR, IPXS, ISTTF, IWDS, IXXY, KSS, KZEF,
     *                  N, NBF, NBFQ, NFR, NSTTF, NXSP
C     .. Array Arguments ..
      DOUBLE PRECISION  BF(NBFQ), PXS(IPXS,NXSP), STTF(ISTTF), W(IDW),
     *                  WDS(IWDS,NXSP), XXY(IXXY,NXSP)
      INTEGER           MR(7), MT(4,NXSP)
C     .. Local Scalars ..
      DOUBLE PRECISION  ZERO
      INTEGER           I, IFAILQ, J, JT, JW, KF1, KF10, KF11, KF12,
     *                  KF13, KF2, KF3, KF4, KF5, KF6, LEA, LER, LEW,
     *                  NB, ND, NDD, NDS, NDV, NEW, NEX, NGW, NMBQ, NMP,
     *                  NNB, NNP, NNQ, NNR, NP, NPAR, NPARQ, NPD, NPS,
     *                  NPX, NQ, NQD, NQQ, NQS, NS, NU, NWD, NXS
C     .. Local Arrays ..
      INTEGER           MPQS(4)
C     .. External Subroutines ..
      EXTERNAL          G13AAF, G13AEU, G13AJY, G13BEK, G13BEL, G13BEX
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Data statements ..
      DATA              ZERO/0.0D0/
C     .. Executable Statements ..
C
C     DERIVE MISCELLANEOUS INTEGER VALUES ASSOCIATED WITH Y ORDERS
C
      CALL G13AJY(MR,NP,ND,NQ,NPS,NDS,NQS,NS,NPD,NDD,NQD,MPQS,NPAR)
C
C     DERIVE COMPONENT LENGTHS OF STATE SET
C
      LEW = NS*NPS + NDD
      LEA = MAX((NS*NQS),NP)
      LER = NQ
      IERR = 11
      NSTTF = 0
      NPARQ = MAX(NPAR,1)
C
C     DERIVE START POINTS OF COMPONENT PARTS OF WORKING ARRAY
C
      KF1 = 1
      KF2 = KF1 + N
      KF3 = KF2 + N
      KF4 = KF3 + NFR
      KF5 = KF4 + NFR
      KF6 = KF5 + IPXS
      KF10 = KF6 + KSS
      KF11 = KF10 + NPARQ
      KF12 = KF11 + NPARQ
      KF13 = KF12 + NPARQ
C
C     START TO BUILD UP E SERIES
C
      DO 20 JT = 1, N
         W(JT) = XXY(JT,NXSP)
   20 CONTINUE
      IF (NXSP.LE.1) GO TO 200
      NXS = NXSP - 1
C
C     PROCESS EACH INPUT SERIES IN TURN
C
      DO 180 I = 1, NXS
         CALL G13BEX(MT,I,NXSP,NNB,NNP,NNQ,NNR,NWD,NGW,NPX)
         NMBQ = N
         NMP = N
         IF (NNR.GT.3 .OR. NNR.LT.2) GO TO 40
         NMBQ = N - NNB - NNQ
         NMP = N - NNP
C
C        TRANSFER (B+Q) VALUES OF X FROM THIS INPUT SERIES TO
C        THE STATE SET
C
   40    JW = KF2 - 1
         DO 60 JT = 1, N
            JW = JW + 1
            W(JW) = XXY(JT,I)
            IF (JT.LE.NMBQ) GO TO 60
            NSTTF = NSTTF + 1
            IF (NSTTF.GT.ISTTF) GO TO 400
            STTF(NSTTF) = W(JW)
   60    CONTINUE
C
C        DERIVE THE Z(T) RESPONSES FOR THIS INPUT SERIES
C
         IF (NPX.LE.0) GO TO 100
         JW = KF5 - 1
         DO 80 J = 1, NPX
            JW = JW + 1
            W(JW) = PXS(J,I)
   80    CONTINUE
  100    JW = KF6 - 1
         DO 120 J = 1, NWD
            JW = JW + 1
            W(JW) = WDS(J,I)
  120    CONTINUE
         CALL G13BEL(W(KF2),N,W(KF5),IPXS,W(KF6)
     *               ,NWD,NNB,NNP,NNQ,NNR,NPX,NEX,W(KF3),W(KF4),NFR)
         JW = KF4 + NPX - 1
C
C        TRANSFER P VALUES OF Z(T) TO THE STATE SET. BUILD UP E IN W(1).
C        REPLACE X,Y WITH Z,E IF KZEF.NE.0
C
         DO 160 JT = 1, N
            JW = JW + 1
            IF (JT.LE.NMP) GO TO 140
            NSTTF = NSTTF + 1
            IF (NSTTF.GT.ISTTF) GO TO 400
            STTF(NSTTF) = W(JW)
  140       W(JT) = W(JT) - W(JW)
            IF (KZEF.EQ.0) GO TO 160
            XXY(JT,I) = W(JW)
            XXY(JT,NXSP) = XXY(JT,NXSP) - W(JW)
  160    CONTINUE
  180 CONTINUE
  200 IFAILQ = 1
C
C     DIFFERENCE THE E SERIES
C
      CALL G13AAF(W,N,ND,NDS,NS,W(KF3),NDV,IFAILQ)
      GO TO 260
  220 IF (C.EQ.ZERO) GO TO 300
C
C     SUBTRACT C FROM THE DIFFERENCED E SERIES TO GIVE THE W SERIES.
C     THE RECONSTITUTION DATA REMAINS UNCHANGED.
C
      JW = KF3 - 1
      DO 240 I = 1, NDV
         JW = JW + 1
         W(JW) = W(JW) - C
  240 CONTINUE
      GO TO 300
  260 IF (LEW.LE.0) GO TO 220
C
C     PUT DIFFERENCED E SERIES AND RECONSTITUTION DATA
C     INTO THE STATE SET
C
      JW = KF3 - 1 + N - LEW
      DO 280 I = 1, LEW
         JW = JW + 1
         NSTTF = NSTTF + 1
         IF (NSTTF.GT.ISTTF) GO TO 400
         STTF(NSTTF) = W(JW)
  280 CONTINUE
      GO TO 220
C
C     ADD BACK FORECASTS TO THE FRONT OF THE W SERIES
C
  300 CALL G13BEK(BF,NBFQ,W(KF3),NFR,W(KF3),NFR,NBF,N,NQQ)
      NEW = NQQ - NDD
      NB = 0
      NU = 1
C
C     DERIVE A AND ALPHA SERIES
C
      CALL G13AEU(0,W(KF3),W(KF4),W(1),NEW,W(KF6),W(KF6+1),W(KF6-1)
     *            ,NU,W(KF10),W(KF11),W(KF12),W(KF13)
     *            ,NPARQ,NP,NQ,NPS,NQS,NS,NB)
      IF (LEA.LE.0) GO TO 340
C
C     TRANSFER ALPHA VALUES TO THE STATE SET
C
      JW = KF4 - 1 + NEW - LEA
      DO 320 I = 1, LEA
         JW = JW + 1
         NSTTF = NSTTF + 1
         IF (NSTTF.GT.ISTTF) GO TO 400
         STTF(NSTTF) = W(JW)
  320 CONTINUE
  340 IF (LER.LE.0) GO TO 380
C
C     TRANSFER A VALUES TO THE STATE SET
C
      JW = NEW - LER
      DO 360 I = 1, LER
         JW = JW + 1
         NSTTF = NSTTF + 1
         IF (NSTTF.GT.ISTTF) GO TO 400
         STTF(NSTTF) = W(JW)
  360 CONTINUE
  380 IERR = 0
  400 RETURN
      END
