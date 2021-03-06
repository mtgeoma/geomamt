      SUBROUTINE G13BGF(STTF,NSTTF,MR,NSER,MT,PARA,NPARA,NNV,XXYN,IXXYN,
     *                  KZEF,RES,WA,IWA,IFAIL)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11A REVISED. IER-453 (JUN 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G13BFF UPDATES A STATE SET WHEN NEW VALUES OF X AND Y
C     ARE ADDED TO A FULLY SPECIFIED TRANSFER FUNCTION TIME
C     SERIES MODEL
C
C
C     BREAK DOWN MR(ORDERS ARRAY FOR Y) INTO COMPONENT PARTS
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13BGF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IWA, IXXYN, KZEF, NNV, NPARA, NSER, NSTTF
C     .. Array Arguments ..
      DOUBLE PRECISION  PARA(NPARA), RES(NNV), STTF(NSTTF), WA(IWA),
     *                  XXYN(IXXYN,NSER)
      INTEGER           MR(7), MT(4,NSER)
C     .. Local Scalars ..
      INTEGER           I, IERROR, KQ, LPARA, LSTTF, ND, NDD, NDS, NGW,
     *                  NNB, NNP, NNQ, NNR, NP, NPAR, NPD, NPS, NPX, NQ,
     *                  NQD, NQS, NS, NWD, NXS
C     .. Local Arrays ..
      INTEGER           MPQS(4)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G13AJY, G13BEX, G13BFZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
      CALL G13AJY(MR,NP,ND,NQ,NPS,NDS,NQS,NS,NPD,NDD,NQD,MPQS,NPAR)
C
C     USE INFORMATION IN MR(ORDERS ARRAY FOR Y) AND MT(ORDERS
C     ARRAY FOR XS) TO CALCULATE NECESSARY LENGTHS OF PARA
C     AND STTF
C
      LSTTF = NS*NPS + NDD + NQ + MAX(NP,NS*NQS)
      LPARA = NPAR + 1
      IF (NSER.LE.1) GO TO 40
      NXS = NSER - 1
      IERROR = 5
      DO 20 I = 1, NXS
         IF (MT(4,I).LT.1 .OR. MT(4,I).GT.3) GO TO 60
         CALL G13BEX(MT,I,NSER,NNB,NNP,NNQ,NNR,NWD,NGW,NPX)
         LSTTF = LSTTF + NNB + NNQ + NNP
         LPARA = LPARA + NWD
   20 CONTINUE
C
C     RETURN WITH RELEVANT IFAIL VALUE IF AN ERROR IS FOUND
C
   40 IERROR = 1
      IF (NSTTF.NE.LSTTF) GO TO 60
      IERROR = 2
      IF (NPARA.NE.LPARA) GO TO 60
      IERROR = 3
      IF (NNV.GT.IXXYN) GO TO 60
      KQ = NNV + 2*NSTTF + MAX(NNV,4*NPAR) + MAX(NNV,NSTTF)
      IERROR = 4
      IF (IWA.LT.KQ) GO TO 60
C
C     CALL AUXILIARY ROUTINE TO DO UPDATE CALCULATIONS
C
      CALL G13BFZ(STTF,NSTTF,MT,NSER,PARA,NPARA,MR,XXYN,IXXYN,NNV,KZEF,
     *            RES,WA,IWA)
      IFAIL = 0
      RETURN
   60 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
