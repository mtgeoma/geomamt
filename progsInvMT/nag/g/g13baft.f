      SUBROUTINE G13BAF(Y,NY,MR,NMR,PAR,NPAR,CY,WA,NWA,B,NB,IFAIL)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G13BAF FILTERS A TIME SERIES BY AN AUTOREGRESSIVE
C     INTEGRATED MOVING AVERAGE MODEL.
C
C     CONTRIBUTORS - G.TUNNICLIFFE WILSON,M. HURLEY (LANC. UNIV.)
C     VALIDATOR    - T. LAMBERT ( NAG CENTRAL OFFICE )
C
C     USES NAG LIBRARY ROUTINES G13BAY, G13BAZ, P01AAF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13BAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CY
      INTEGER           IFAIL, NB, NMR, NPAR, NWA, NY
C     .. Array Arguments ..
      DOUBLE PRECISION  B(NB), PAR(NPAR), WA(NWA), Y(NY)
      INTEGER           MR(NMR)
C     .. Local Scalars ..
      DOUBLE PRECISION  EF, W0
      INTEGER           I, IB, IB1, IB2, IDD, IDYD, IERR, IERROR, IPD,
     *                  IPYD, IQD, IQYD, M, MY, N
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G13BAY, G13BAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
      IERROR = 1
C
C     CHECK LENGTH OF ORDERS ARRAY
C
      IF (NMR.NE.7 .AND. NMR.NE.14) GO TO 140
C
C     CHECK ORDERS ARRAY
C
      IERROR = 2
      DO 20 I = 1, NMR
         IF (MR(I).LT.0) GO TO 140
   20 CONTINUE
      M = MR(1) + MR(3) + MR(4) + MR(6)
      IF (M.LE.0) GO TO 140
      IF (MR(7).EQ.1) GO TO 140
      N = MR(4) + MR(5) + MR(6)
      IF (MR(7).EQ.0 .AND. N.GT.0) GO TO 140
      IF (MR(7).NE.0 .AND. N.EQ.0) GO TO 140
      IF (NMR.EQ.7) GO TO 40
      IF (MR(14).EQ.1) GO TO 140
      N = MR(11) + MR(12) + MR(13)
      IF (MR(14).EQ.0 .AND. N.GT.0) GO TO 140
      IF (MR(14).NE.0 .AND. N.EQ.0) GO TO 140
   40 CONTINUE
      IQD = MR(3) + MR(6)*MR(7)
      IDD = MR(2) + MR(5)*MR(7)
      IPD = MR(1) + MR(4)*MR(7)
      IF (NMR.EQ.7) GO TO 60
      IQYD = MR(10) + MR(13)*MR(14)
      IDYD = MR(9) + MR(12)*MR(14)
      IPYD = MR(8) + MR(11)*MR(14)
   60 CONTINUE
C
C     CHECK NUMBER OF PARAMETERS
C
      IERROR = 3
      MY = 0
      IF (NMR.EQ.14) MY = MR(8) + MR(10) + MR(11) + MR(13)
      IF (NPAR.NE.M+MY) GO TO 140
C
C     CHECK LENGTH OF SERIES
C
      IERROR = 4
      N = 0
      IF (NMR.EQ.14) N = IQYD
      N = MAX(N+1,NPAR)
      IF (NY.LT.N) GO TO 140
      IB1 = M
      IB2 = MY
C
C     CHECK SUFFICIENT WORKSPACE
C
      IERROR = 5
      M = 1
      IF (NMR.EQ.7) GO TO 80
      M = IQD + IPYD + IDYD
      M = M*(M+2)
   80 CONTINUE
      IF (NWA.LT.M) GO TO 140
C
C     CHECK DIMENSION OF B
C
      IERROR = 6
      M = NY
      IF (NMR.EQ.14) M = M + MAX(IQD,IDD+IPD)
      IF (NB.LT.M) GO TO 140
C
C     CHECK PARAMETERS VALID
C
      IERROR = 7
      EF = 1000.0D0
      DO 100 I = 1, NPAR
         B(I) = PAR(I)
  100 CONTINUE
      CALL G13BAY(MR,B,IB1,EF,IERR)
      IF (IERR.NE.0) GO TO 140
      IF (NMR.EQ.7) GO TO 120
      IERROR = 8
      IB1 = IB1 + 1
      IF (IB2.EQ.0) GO TO 120
      CALL G13BAY(MR(8),B(IB1),IB2,EF,IERR)
      IF (IERR.NE.0) GO TO 140
  120 CONTINUE
C     IF THIS POINT IS REACHED PARAMETER
C     VALUES HAVE PASSED INITIAL CHECKS
      IERROR = 0
      W0 = 1.0D0
      IB = -1
      CALL G13BAZ(Y,NY,MR,NMR,IB,W0,PAR,NPAR,CY,WA,NWA,B,NB,IERROR)
      IF (IERROR.NE.0) GO TO 140
C     SUCCESSFUL EXIT
      IFAIL = 0
      RETURN
C     UNSUCCESSFUL EXIT
  140 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
