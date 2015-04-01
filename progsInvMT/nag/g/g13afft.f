      SUBROUTINE G13AFF(MR,PAR,NPAR,C,KFC,X,NX,S,NDF,SD,NPPC,CM,ICM,ST,
     *                  NST,KPIV,NIT,ITC,ISF,RES,IRES,NRES,IFAIL)
C     MARK 9 RELEASE. NAG COPYRIGHT 1981.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G13AFF IS AN EASY TO USE VERSION OF G13AEF.
C     IT FITS A SEASONAL ARIMA MODEL TO AN OBSERVED TIME SERIES,
C     USING
C     A NON-LINEAR LEAST SQUARES TECHNIQUE INCORPORATING
C     BACKFORECASTING
C
C     USES NAG LIBRARY ROUTINES G13AEF AND P01AAF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13AFF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  C, S
      INTEGER           ICM, IFAIL, IRES, ITC, KFC, KPIV, NDF, NIT,
     *                  NPAR, NPPC, NRES, NST, NX
C     .. Array Arguments ..
      DOUBLE PRECISION  CM(ICM,NPPC), PAR(NPAR), RES(IRES), SD(NPPC),
     *                  ST(NX), X(NX)
      INTEGER           ISF(4), MR(7)
C     .. Local Scalars ..
      INTEGER           I, IERROR, IEX, IGH, IH, IJ, IQ, IRESQ, IWA, J,
     *                  JQ, KAL, KEX, KEXR, KG, KH, KHC, KSD, KWA, KZSP,
     *                  ND, NDS, NEX, NGH, NP, NPD, NPS, NQ, NQD, NQS,
     *                  NS
C     .. Local Arrays ..
      INTEGER           ICOUNT(6)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G13AEF, G13AFZ
C     .. Executable Statements ..
      NP = MR(1)
      ND = MR(2)
      NQ = MR(3)
      NPS = MR(4)
      NDS = MR(5)
      NQS = MR(6)
      NS = MR(7)
      NPD = NP + NPS*NS
      NQD = NQ + NQS*NS
C
C     CHECK ORDERS VECTOR, NPAR, AND NPPC FOR CONSISTENCY
C
      IERROR = 1
      IF (KFC.NE.0 .AND. KFC.NE.1) GO TO 160
      DO 20 I = 1, 7
         IF (MR(I).LT.0) GO TO 160
   20 CONTINUE
      IF (NP+NQ+NPS+NQS.LT.1) GO TO 160
      IF (NP+NQ+NPS+NQS.NE.NPAR) GO TO 160
      IF (NPAR+KFC.NE.NPPC) GO TO 160
      IF (NS.EQ.1) GO TO 160
      IF (NS.EQ.0 .AND. NPS+NDS+NQS.NE.0) GO TO 160
      IF (NS.NE.0 .AND. NPS+NDS+NQS.EQ.0) GO TO 160
C
C     CHECK THAT ICM IS LARGE ENOUGH
C
      IERROR = 6
      IF (ICM.LT.NPPC) GO TO 160
C
C     CALCULATE WORKSPACE REQUIREMENTS AND CHECK SIZE OF RES ARRAY
C
      NGH = NQD + NPPC
      IRESQ = 15*NQD + 11*NX + 13*NPPC + 8*NPD + 12 + 2*NGH*NGH
      IERROR = 5
      IF (IRES.LT.IRESQ) GO TO 160
C
C     G13AEF CAN NOW BE CALLED SAFELY AND SHOULD DETECT
C     ANY REMAINING ERRORS IN THE G13AFF PARAMETERS
C
C     FIRST FORM START POINTS OF COMPONENT ARRAYS IN RES
C
      NEX = NX + NQD
      NQ = (NGH+1)*NGH
      KEX = 1
      KEXR = NEX + 1
      KAL = KEXR + NEX
      KG = KAL + NEX
      KSD = KG + NGH
      KH = KSD + NGH
      KHC = KH + NQ
      KZSP = KHC + NQ
      KWA = KZSP + 4
      IWA = IRES - KWA
      IGH = NGH
      IEX = NEX
      IH = IGH + 1
C
C     CARRY OUT ESTIMATION
C
      IERROR = 1
      CALL G13AEF(MR,PAR,NPAR,C,KFC,X,NX,ICOUNT,RES(KEX),RES(KEXR)
     *            ,RES(KAL),IEX,S,RES(KG),IGH,RES(KSD),RES(KH)
     *            ,IH,ST,NX,NST,G13AFZ,KPIV,NIT,ITC,RES(KZSP)
     *            ,0,ISF,RES(KWA),IWA,RES(KHC),IERROR)
      NDF = ICOUNT(5)
      IF (IERROR.EQ.0) GO TO 40
C
C     G13AEF FAILS - CHECK FOR UNEXPECTED FAILURES
C
      IF (IERROR.EQ.5 .OR. IERROR.EQ.6) IERROR = 11
      GO TO 160
C
C     G13AEF EXITS SUCCESSFULLY
C
C     COPY STANDARD DEVIATIONS, CORRELATIONS, AND RESIDUALS
C     FROM WORK ARRAY INTO SD, CM, AND RES
C
   40 J = NQD + KSD - 1
      DO 60 I = 1, NPPC
         J = J + 1
         SD(I) = RES(J)
   60 CONTINUE
      DO 120 J = 1, NPPC
         JQ = NQD + J
         DO 100 I = J, NPPC
            IF (J.NE.I) GO TO 80
            CM(I,I) = 1.0D0
            GO TO 100
   80       IQ = NQD + I
            IJ = KH - 1 + IQ + IH*(JQ-1)
            CM(I,J) = RES(IJ)
            CM(J,I) = RES(IJ)
  100    CONTINUE
  120 CONTINUE
      J = NEX + NQD
      NRES = ICOUNT(2)
      DO 140 I = 1, NRES
         J = J + 1
         RES(I) = RES(J)
  140 CONTINUE
      IFAIL = 0
      RETURN
  160 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
