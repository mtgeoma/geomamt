      SUBROUTINE D03PUF(ULEFT,URIGHT,GAMMA,FLUX,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     This routine implements the Roe Approximate Riemann Solver for
C     the Euler equations, primarily for use with D03PFF/D03PLF/D03PSF.
C     New version, 10/8/94
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D03PUF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  GAMMA
      INTEGER           IFAIL
C     .. Array Arguments ..
      DOUBLE PRECISION  FLUX(3), ULEFT(3), URIGHT(3)
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, A2, A3, AL1, AL2, AL3, C, D1, D2, D3, G1,
     *                  HA, HL, HR, HUA, L1, L2, L3, PL, PR, R, RR, UA,
     *                  UL, UR
      INTEGER           I, IFAIL1, J
      CHARACTER*200     ERRMSG
C     .. Local Arrays ..
      DOUBLE PRECISION  E(3,3), FL(3), FR(3)
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02NNQ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      IFAIL1 = 0
      IF (GAMMA.LE.0.0D0) GO TO 60
C
C     Deltas ..
      D1 = URIGHT(1) - ULEFT(1)
      D2 = URIGHT(2) - ULEFT(2)
      D3 = URIGHT(3) - ULEFT(3)
C
C     Check for negative densities ..
      IF (ULEFT(1).LT.0.0D0 .OR. URIGHT(1).LT.0.0D0) GO TO 80
C
C     Factors for averages ..
      R = SQRT(URIGHT(1)/ULEFT(1))
      RR = 1.0D0/(1.0D0+R)
C
C     Left and right values of u, p and h ..
      UL = ULEFT(2)/ULEFT(1)
      UR = URIGHT(2)/URIGHT(1)
      G1 = GAMMA - 1.0D0
      PL = G1*(ULEFT(3)-0.5D0*ULEFT(2)**2/ULEFT(1))
      PR = G1*(URIGHT(3)-0.5D0*URIGHT(2)**2/URIGHT(1))
C     Check for negative pressures ..
      IF (PL.LT.0.0D0 .OR. PR.LT.0.0D0) GO TO 100
      HL = (ULEFT(3)+PL)/ULEFT(1)
      HR = (URIGHT(3)+PR)/URIGHT(1)
C
C     Average values of u and h ..
      UA = RR*(UL+R*UR)
      HA = RR*(HL+R*HR)
      HUA = HA - 0.5D0*UA**2
C     Sound speed ..
      C = SQRT(G1*HUA)
C
C     Eigenvalues ..
      L1 = UA - C
      L2 = UA
      L3 = UA + C
C     Eigenvectors ..
C     (E(I,J) contains Ith component of Jth e-vector)
      DO 20 J = 1, 3
         E(1,J) = 1.0D0
   20 CONTINUE
      E(2,1) = UA - C
      E(2,2) = UA
      E(2,3) = UA + C
      E(3,1) = HA - UA*C
      E(3,2) = 0.5D0*UA**2
      E(3,3) = HA + UA*C
C
C     Alphas ..
      A2 = -(D3+D1*(UA**2-HA)-UA*D2)/HUA
      A3 = 0.5D0*(D2+D1*(C-UA)-A2*C)/C
      A1 = D1 - A2 - A3
C
C     Left and right fluxes ..
      FL(1) = ULEFT(2)
      FL(2) = PL + ULEFT(1)*UL**2
      FL(3) = UL*(ULEFT(3)+PL)
      FR(1) = URIGHT(2)
      FR(2) = PR + URIGHT(1)*UR**2
      FR(3) = UR*(URIGHT(3)+PR)
C     Roe fluxes ..
      AL1 = A1*ABS(L1)
      AL2 = A2*ABS(L2)
      AL3 = A3*ABS(L3)
      DO 40 I = 1, 3
         FLUX(I) = 0.5D0*(FL(I)+FR(I)) - 0.5D0*(AL1*E(I,1)+AL2*E(I,2)
     *             +AL3*E(I,3))
   40 CONTINUE
      GO TO 120
C
   60 CONTINUE
      ERRMSG =
     *' Routine entered with GAMMA (=R1) less than or equal     to 0.0.
     *'
      CALL D02NNQ(ERRMSG,1,0,0,0,1,GAMMA,0.0D0)
      IFAIL1 = 1
      GO TO 120
C
   80 CONTINUE
      ERRMSG =
     *' Routine entered with left density value i.e. ULEFT(1)   (=R1) an
     *d/or right density value i.e. URIGHT(1) (=R2) less than   0.0. '
      CALL D02NNQ(ERRMSG,1,0,0,0,2,ULEFT(1),URIGHT(1))
      IFAIL1 = 2
      GO TO 120
C
  100 CONTINUE
      ERRMSG =
     *' Routine entered with left pressure value (=R1) and/or   right pr
     *essure value (=R2) less than 0.0. '
      CALL D02NNQ(ERRMSG,1,0,0,0,2,PL,PR)
      IFAIL1 = 2
C
  120 CONTINUE
      IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,0,REC)
      RETURN
      END
