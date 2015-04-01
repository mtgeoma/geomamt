      SUBROUTINE D03PDM(NPDE,NPTS,U,X,OMEGA,DU,XBK,NEL,NPTL,XC,CCR,XBH,
     *                  IBK,DUTEM,V,NV,AUXINI,DINPDF,DINPJF)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C  ---------------------------------------------------------------------
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C       Fortran functions used:  sin cos .
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  ---------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IBK, NEL, NPDE, NPTL, NPTS, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  CCR(NPTL), DU(NPTL,NPTL), DUTEM(NPTL,NPTL),
     *                  OMEGA(NPTL,NPTL), U(NPDE,NPTS), V(*), X(NPTS),
     *                  XBH(IBK), XBK(IBK), XC(NPTL)
C     .. Subroutine Arguments ..
      EXTERNAL          AUXINI, DINPDF, DINPJF
C     .. Arrays in Common ..
      DOUBLE PRECISION  CCRULE(50)
C     .. Local Scalars ..
      DOUBLE PRECISION  H1, H2, PI, SINT, SUM, TEMP
      INTEGER           I, IJ, ITEM, J, K, NM1, NT, NTP1
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE, SIN
C     .. Common blocks ..
      COMMON            /PD03PD/CCRULE
C     .. Save statement ..
      SAVE              /PD03PD/
C     .. Executable Statements ..
C
C ... Form  constants for wkspace initialisation ...
C
      NM1 = NPTL - 1
      PI = X01AAF(0.0D0)
C
C ... Formation of grid and initial values of U ...
C
      DO 40 I = 1, NEL
         H1 = XBH(I+1) - XBH(I)
         H2 = XBH(I+1) + XBH(I)
         XBK(I) = XBH(I)
         DO 20 J = 1, NPTL
            IJ = (I-1)*NM1 + J
            IF (I.EQ.1) XC(J) = COS(PI*DBLE(J-NPTL)/NM1)
            X(IJ) = (XC(J)*H1+H2)*0.5D0
            IF (J.EQ.1) X(IJ) = XBH(I)
            IF (J.EQ.NPTL) X(IJ) = XBH(I+1)
   20    CONTINUE
   40 CONTINUE
      XBK(IBK) = XBH(IBK)
      XC(1) = -1.0D0
      XC(NPTL) = 1.0D0
C
C ... Form the matrix OMEGA ...
C
      DO 80 J = 1, NPTL
         DO 60 I = 1, NPTL
            OMEGA(I,J) = 2.D0*COS(PI*(I-1)*(NPTL-J)/NM1)/NM1
   60    CONTINUE
   80 CONTINUE
C
C ... Modify edges of OMEGA and form edges ...
C ... of intermediate DU matrix            ...
C
      ITEM = 1
      DO 100 I = 1, NPTL
         OMEGA(I,1) = OMEGA(I,1)*0.5D0
         OMEGA(1,I) = OMEGA(1,I)*0.5D0
         OMEGA(NPTL,I) = OMEGA(NPTL,I)*0.5D0
         OMEGA(I,NPTL) = OMEGA(I,NPTL)*0.5D0
         DUTEM(I,1) = 0.0D0
         DUTEM(1,I) = -DBLE((I-1)**2*ITEM)
         DUTEM(NPTL,I) = DBLE((I-1)**2)
         ITEM = -ITEM
  100 CONTINUE
C
C ... Finish forming rest of intermediate DU matrix that is ...
C ... held in DUTEM                                         ...
C
      IF (NPTL.GT.2) THEN
         DO 140 I = 2, NM1
            TEMP = PI*(I-NPTL)/NM1
            SINT = SIN(TEMP)
            DO 120 J = 2, NM1
               DUTEM(I,J) = SIN(TEMP*(J-1))/SINT*(J-1)
  120       CONTINUE
            DUTEM(I,NPTL) = 0.0D0
  140    CONTINUE
      END IF
C
C ... Form full DU by matrix multiplication ...
C
      DO 200 I = 1, NPTL
         DO 180 J = 1, NPTL
            DU(I,J) = 0.0D0
            DO 160 K = 1, NPTL
               DU(I,J) = DU(I,J) + DUTEM(I,K)*OMEGA(K,J)
  160       CONTINUE
  180    CONTINUE
  200 CONTINUE
C
C ... Calculate the coeffs of the Clenshaw-Curtis rule ...
C
      NT = NM1/2
      IF ((2*NT).NE.NM1) NT = (NM1-1)/2
      NTP1 = NT + 1
      SUM = 0.0D0
      DO 240 I = 1, NPTL
         TEMP = 0.5D0
         CCR(I) = 0.0D0
         DO 220 K = 1, NTP1
            IF (K.EQ.NTP1 .AND. ((2*NT).EQ.NM1)) TEMP = 0.5D0
            CCR(I) = CCR(I) + COS(2.0D0*(I-1)*(K-1)*PI/NM1)
     *               *TEMP/(4.0D0*(K-1)**2-1.0D0)
            TEMP = 1.0D0
  220    CONTINUE
         IF (I.EQ.1 .OR. I.EQ.NPTL) TEMP = 0.5D0
         CCR(I) = CCR(I)*(-4.0D0)*TEMP/NM1
         SUM = SUM + CCR(I)
  240 CONTINUE
      DO 260 I = 1, NPTL
         CCRULE(I) = CCR(I)
  260 CONTINUE
      DO 280 I = 2, NM1
         CCR(I) = CCR(I)/CCR(1)
  280 CONTINUE
C
C ... Find the initial values of the ODE and PDE components ...
C
      CALL AUXINI(DINPDF,DINPJF,NPDE,NPTS,X,U,NV,V)
      RETURN
C
      END
