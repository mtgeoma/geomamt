      SUBROUTINE G13AEV(NP,NQ,NPS,NQS,NS,AQ,NAQ,BQ,NBQ,KSCH,AL,NA,NB,S,
     *                  G,H,IZ,NGH)
C     MARK 9 RELEASE. NAG COPYRIGHT 1981.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G13AEV CALCULATES S, G AND H, WHICH ARE EQUIVALENT TO THE
C     RESIDUAL SUM OF SQUARES SYY, AND TO SYX AND SXX RESPECTIVELY
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  S
      INTEGER           IZ, KSCH, NA, NAQ, NB, NBQ, NGH, NP, NPS, NQ,
     *                  NQS, NS
C     .. Array Arguments ..
      DOUBLE PRECISION  AL(NA), AQ(NAQ), BQ(NBQ), G(NGH), H(IZ,NGH)
C     .. Local Scalars ..
      DOUBLE PRECISION  ZERO
      INTEGER           I, IH, IV, J, JH, JHL, JV, KHA, KHAQ, KHB, KHBQ,
     *                  KHS, KVA, KVAQ, KVB, KVBQ, KVS, NAL, NGHP, NHL,
     *                  NID, NPA, NPD, NQD, NVL
C     .. Local Arrays ..
      INTEGER           MAK(7), MAL(7), MAS(7), MBK(7), MBL(7), MNC(7)
C     .. Data statements ..
      DATA              MAK(1)/0/, MAK(2)/1/, MAK(3)/1/, MAK(4)/1/,
     *                  MAK(7)/0/
      DATA              MAL(1)/0/, MAL(2)/-1/, MAL(3)/0/, MAL(4)/0/,
     *                  MAL(5)/0/, MAL(6)/0/
      DATA              MBK(1)/0/, MBK(2)/1/, MBK(3)/-1/, MBK(4)/1/,
     *                  MBK(7)/0/
      DATA              MBL(1)/0/, MBL(2)/-1/, MBL(4)/0/, MBL(6)/0/,
     *                  MBL(7)/0/
      DATA              MAS(1)/1/, MAS(2)/-1/, MAS(3)/1/, MAS(4)/-1/,
     *                  MAS(5)/1/
      DATA              MAS(6)/-1/, MAS(7)/1/, MAL(7)/0/, MNC(1)/1/,
     *                  MNC(7)/1/
      DATA              ZERO/0.0D0/
C     .. Executable Statements ..
      DO 40 J = 1, NGH
         G(J) = ZERO
         DO 20 I = 1, IZ
            H(I,J) = ZERO
   20    CONTINUE
   40 CONTINUE
      NPD = NP + NPS*NS
      NQD = NQ + NQS*NS
C
C     ASSIGN INITIAL VALUES NOT COVERED BY DATA STATEMENTS
C
      MAK(5) = NS
      MAK(6) = NS
      MBK(5) = -NS
      MBK(6) = NS
      MBL(3) = NP
      MBL(5) = NPS*NS
      MNC(2) = NQD
      MNC(3) = NP
      MNC(4) = NQ
      MNC(5) = NPS
      MNC(6) = NQS
C
C     NID HOLDS THE NUMBER OF A AND B SETS IN AQ AND BQ WHICH ARE
C     NEEDED IN THE CALCULATION OF S, G AND H
C
      NID = KSCH
      IF (NID.LE.2) GO TO 60
      NID = NID + 3
   60 NPA = 0
      NGHP = NGH + 1
C
C     LOOP OVER EACH A AND B SET, AS THOUGH THESE RELATED TO SETS
C     OF ROWS IN A TWO DIMENSIONAL TABLE HOLDING NGHP ROWS AND
C     COLUMNS, WHERE NGHP IS THE NUMBER OF ESTIMATED PARAMETERS
C     PLUS 1
C
      DO 380 IV = 1, NID
         IF (MNC(IV).LE.0) GO TO 380
         NVL = MNC(IV)
C
C        THE MAS VALUES INDICATE THE PRESENCE OF A POSITIVE OR NEGATIVE
C        MULTIPLYING FACTOR. KVS HOLDS THE FACTOR FOR THIS SET OF ROWS
C
         KVS = MAS(IV)
C
C        LOOP OVER THE NUMBER OF ESTIMATED PARAMETERS IN THIS SET, AS
C        THOUGH THESE RELATED TO THE SEPARATE ROWS
C
         DO 360 JV = 1, NVL
C
C           KVA AND KVB ARE USED TO DETERMINE THE START POINT OF
C           THE LOOP USED IN THE CALCULATION OF EACH AL VALUE
C
            KVA = MAK(IV)*JV + MAL(IV)
            KVB = MBK(IV)*JV + MBL(IV)
C
C           ZEROISE THE AL VALUES RELATING TO THIS ROW
C
            DO 80 I = 1, NGHP
               AL(I) = ZERO
   80       CONTINUE
C
C           NAL WILL HOLD THE COUNT OF AL VALUES ALONG THIS ROW, AND NPA
C           WILL HOLD THE ROW COUNT
C
            NAL = NPA
            NPA = NPA + 1
C
C           LOOP OVER A AND B SETS AS THOUGH THESE RELATED TO SETS OF
C           COLUMNS TO THE RIGHT OF THE MAIN DIAGONAL
C
            DO 260 IH = IV, NID
               IF (MNC(IH).LE.0) GO TO 260
               NHL = MNC(IH)
C
C              KHS HOLDS THE MULTIPLYING FACTOR FOR THIS SET OF COLUMNS
C
               KHS = MAS(IH)
               JHL = JV
               IF (IH.EQ.IV) GO TO 100
               JHL = 1
C
C              LOOP OVER CERTAIN OF THE ESTIMATED PARAMETERS IN
C              THIS SET, AS THOUGH THESE RELATED TO THE SEPARATE SETS
C              OF COLUMNS TO THE RIGHT OF THE MAIN DIAGONAL
C
  100          DO 240 JH = JHL, NHL
C
C                 KHA AND KHB ARE USED TO DETERMINE THE START POINT
C                 OF THE LOOP USED IN THE CALCULATION OF EACH AL VALUE
C
                  KHA = MAK(IH)*JH + MAL(IH)
                  KHB = MBK(IH)*JH + MBL(IH)
                  NAL = NAL + 1
C
C                 LOOP OVER THE VALUES OF AQ RELATING TO THIS ROW
C                 AND COLUMN
C
                  DO 160 I = 1, NA
                     KVAQ = I - KVA
                     IF (KVAQ.LE.0) GO TO 160
                     KVAQ = KVAQ + NA*(IV-1)
                     KHAQ = I - KHA
                     IF (KHAQ.LE.0) GO TO 160
                     KHAQ = KHAQ + NA*(IH-1)
C
C                    COMPARE THE MULTIPLYING FACTORS AND ADJUST
C                    THE CALCULATIONS OF AL
C
                     IF (KVS-KHS) 120, 140, 120
  120                AL(NAL) = AL(NAL) - AQ(KVAQ)*AQ(KHAQ)
                     GO TO 160
  140                AL(NAL) = AL(NAL) + AQ(KVAQ)*AQ(KHAQ)
  160             CONTINUE
                  IF (NPD.LE.0) GO TO 240
C
C                 LOOP OVER THE VALUES OF BQ RELATING TO THIS ROW
C                 AND COLUMN
C
                  DO 220 I = 1, NPD
                     KVBQ = I - KVB
                     IF (KVBQ.LE.0) GO TO 220
                     KVBQ = KVBQ + NB*(IV-1)
                     KHBQ = I - KHB
                     IF (KHBQ.LE.0) GO TO 220
                     KHBQ = KHBQ + NB*(IH-1)
C
C                    COMPARE THE MULTIPLYING FACTORS AND COMPLETE
C                    THE CALCULATION OF AL
C
                     IF (KVS-KHS) 180, 200, 180
  180                AL(NAL) = AL(NAL) + BQ(KVBQ)*BQ(KHBQ)
                     GO TO 220
  200                AL(NAL) = AL(NAL) - BQ(KVBQ)*BQ(KHBQ)
  220             CONTINUE
  240          CONTINUE
  260       CONTINUE
C
C           S IS GIVEN THE AL VALUE FOR THE FIRST ROW AND COLUMN
C
            IF (IV-1) 320, 280, 320
  280       S = AL(1)
            IF (KSCH.EQ.1) GO TO 400
C
C           THE REMAINING TERMS IN THE FIRST ROW GIVE G
C
            DO 300 I = 2, NAL
               G(I-1) = AL(I)
  300       CONTINUE
            GO TO 360
C
C           AFTER EXCLUDING THE FIRST ROW AND COLUMN, THE VALUES OF AL
C           ARE PUT IN H
C
  320       DO 340 I = NPA, NAL
               H(NPA-1,I-1) = AL(I)
  340       CONTINUE
  360    CONTINUE
  380 CONTINUE
  400 RETURN
      END
