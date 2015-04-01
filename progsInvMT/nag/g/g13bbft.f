      SUBROUTINE G13BBF(Y,NY,MR,NMR,PAR,NPAR,CY,WA,IWA,B,NB,IFAIL)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C        ----------------------------------------------------
C
C        NAG LIBRARY ROUTINE G13BBF FILTERS A TIME SERIES
C        BY AN TF MODEL
C
C        ORIGIN OF SOFTWARE - LANCASTER UNIVERSITY
C        DATE OF INCEPTION - FEBRUARY 1981
C        DATE OF COMPLETION - FEBRUARY 1981
C
C        ALGORITHM DUE TO
C        G.E.P. BOX AND G.M. JENKINS
C        TIME SERIES ANALYSIS  FORECASTING AND CONTROL
C        HOLDEN-DAY
C        BACKFORECASTING FEATURE DUE TO G. TUNICLIFFE-WILSON
C
C        M.A.H.(PROGRAMMER)
C
C        -----------------------------------------------------
C
C        G13BBF CHECKS THE USER SUPPLIED PARAMETERS AND MAKES A
C        FUDGED CALL TO AUXILIARY G13BAZ WHICH CARRIES OUT THE
C        CALCULATIONS.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13BBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CY
      INTEGER           IFAIL, IWA, NB, NMR, NPAR, NY
C     .. Array Arguments ..
      DOUBLE PRECISION  B(NB), PAR(NPAR), WA(IWA), Y(NY)
      INTEGER           MR(NMR)
C     .. Local Scalars ..
      DOUBLE PRECISION  EF
      INTEGER           I, IB, IB1, IB2, IDYD, IERR, IERROR, IPYD, IQYD,
     *                  ISTART, K, L, M, MY, N, NMRD, NPARD, NWAD
C     .. Local Arrays ..
      INTEGER           MRD(14)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G13BAY, G13BAZ, G13BBY, G13BBZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
      IERROR = 0
      IF (IFAIL.NE.0 .AND. IFAIL.NE.1) IFAIL = 0
C        CHECK LENGTH OF ORDERS ARRAY
      IF (NMR.NE.3 .AND. NMR.NE.10) GO TO 220
C        CHECK ORDERS ARRAY
      DO 20 I = 1, NMR
         IF (MR(I).LT.0) GO TO 220
   20 CONTINUE
      IF (NMR.EQ.3) GO TO 40
      IF (MR(10).EQ.1) GO TO 220
      N = MR(7) + MR(8) + MR(9)
      IF (MR(10).EQ.0 .AND. N.GT.0) GO TO 220
      IF (MR(10).NE.0 .AND. N.EQ.0) GO TO 220
   40 CONTINUE
      IF (NMR.EQ.3) GO TO 60
      IQYD = MR(6) + MR(9)*MR(10)
      IDYD = MR(5) + MR(8)*MR(10)
      IPYD = MR(4) + MR(7)*MR(10)
   60 CONTINUE
C        CHECK NUMBER OF PARAMETERS
      M = MR(2) + MR(3) + 1
      MY = 0
      IF (NMR.EQ.10) MY = MR(4) + MR(6) + MR(7) + MR(9)
      IF (NPAR.NE.M+MY) GO TO 220
C        CHECK LENGTH OF SERIES
      N = 0
      IF (NMR.EQ.10) N = IQYD
      N = MAX(N+1,NPAR)
      IF (NY.LT.N) GO TO 260
      IB1 = M
      IB2 = MY
C        CHECK SUFFICIENT WORKSPACE
      M = MR(1) + MR(2) + MR(3) + 1
      L = 0
      IF (NMR.EQ.3) GO TO 80
      N = MR(3) + IPYD + IDYD
      L = N*(N+2)
      M = L + M
   80 M = M + MY
      IF (IWA.LT.M) GO TO 220
C        CHECK DIMENSION OF B
      M = NY
      IF (NMR.EQ.10) M = M + MAX(MR(3),MR(1)+MR(2))
      IF (NB.LT.M) GO TO 220
C        CHECK PARAMETERS VALID
      EF = 1000.0D0
      DO 100 I = 1, NPAR
         B(I) = PAR(I)
  100 CONTINUE
      CALL G13BBZ(MR,B,IB1,EF,IERR)
      IF (IERR.NE.0) GO TO 240
      IF (NMR.EQ.3) GO TO 120
      IF (IB2.EQ.0) GO TO 120
      IB1 = IB1 + 1
      CALL G13BAY(MR(4),B(IB1),IB2,EF,IERR)
      IF (IERR.NE.0) GO TO 240
  120 CONTINUE
C        IF THIS POINT IS REACHED PARAMETER
C        VALUES HAVE PASSED INITIAL CHECKS
C
C        FUDGE PARAMETERS FOR CALL TO G13BAZ
C
      IF (NMR.EQ.3) NMRD = 7
      IF (NMR.EQ.10) NMRD = 14
      MRD(1) = MR(1) + MR(2)
      MRD(2) = 0
      MRD(3) = MR(3)
      MRD(4) = 0
      MRD(5) = 0
      MRD(6) = 0
      MRD(7) = 0
      IF (NMR.EQ.3) GO TO 160
      DO 140 I = 8, 14
         MRD(I) = MR(I-4)
  140 CONTINUE
  160 CONTINUE
      N = L
      ISTART = L
      K = MR(1)
      IF (K.EQ.0) GO TO 200
      DO 180 I = 1, K
         N = N + 1
         WA(N) = 0.0D0
  180 CONTINUE
  200 CONTINUE
      K = MR(2)
      M = 1
      CALL G13BBY(WA,IWA,PAR,NPAR,K,N,M)
      K = MR(3)
      M = MR(2) + 1
      CALL G13BBY(WA,IWA,PAR,NPAR,K,N,M)
      K = IB2
      M = MR(2) + MR(3) + 1
      CALL G13BBY(WA,IWA,PAR,NPAR,K,N,M)
      N = N + 1
      WA(N) = PAR(1)
      NWAD = L
      IF (NWAD.EQ.0) NWAD = 1
      L = L + 1
      NPARD = N - 1 - ISTART
      IB = MR(1)
      CALL G13BAZ(Y,NY,MRD,NMRD,IB,WA(N),WA(L)
     *            ,NPARD,CY,WA,NWAD,B,NB,IERROR)
      IF (IERROR.NE.0) GO TO 280
C        SUCCESSFUL EXIT
      IFAIL = 0
      RETURN
C        UNSUCCESSFUL EXIT
C        ERROR IN USER SUPPLIED PARAMETER
  220 IERROR = 1
      GO TO 280
C        INVALID MODEL(S)
  240 IERROR = 2
      GO TO 280
C        SERIES TOO SHORT
  260 IERROR = 3
  280 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
