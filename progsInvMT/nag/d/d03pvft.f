      SUBROUTINE D03PVF(ULEFT,URIGHT,GAMMA,PATH,FLUX,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     This routine implements the Osher Approximate Riemann Solver
C     for the Euler equations, primarily for use with D03PFF/D03PLF/D03P
C     The parameter PATH on input specifies which variant of the scheme
C     to be used: the P-variant or the O-variant; on output it is unchan
C     unless a sonic point is detected, in which case it will be set to
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D03PVF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  GAMMA
      INTEGER           IFAIL
      CHARACTER         PATH
C     .. Array Arguments ..
      DOUBLE PRECISION  FLUX(3), ULEFT(3), URIGHT(3)
C     .. Local Scalars ..
      DOUBLE PRECISION  AL, C0, C1, C13, C23, CS1, CS2, G0, G1, G2, G3,
     *                  P0, P1, PHI0, PHI1, QSP1, QSP2, R0, R1, U0, U1,
     *                  U13, U23, US1, US2, Z0, Z1
      INTEGER           I, IFAIL1
      CHARACTER         PTH
      CHARACTER*200     ERRMSG
C     .. Local Arrays ..
      DOUBLE PRECISION  D(5), F(3,6), W(3,6)
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02NNQ, E04UDU
C     .. Intrinsic Functions ..
      INTRINSIC         EXP, LOG, SIGN, SQRT
C     .. Executable Statements ..
      IFAIL1 = 0
      PTH = PATH(1:1)
      CALL E04UDU(PTH)
      IF (PTH.NE.'O' .AND. PTH.NE.'P') GO TO 180
      IF (GAMMA.LE.0.0D0) GO TO 200
      G0 = GAMMA - 1.0D0
      G1 = 1.0D0/G0
      G2 = 1/(2.0D0*GAMMA)
      G3 = G0/(GAMMA+1.0D0)
C     Calculate left and right values of primitive variables ..
      R0 = ULEFT(1)
      R1 = URIGHT(1)
C     Check for negative densities ..
      IF (R0.LT.0.0D0 .OR. R1.LT.0.0D0) GO TO 220
      U0 = ULEFT(2)/R0
      U1 = URIGHT(2)/R1
      P0 = G0*(ULEFT(3)-0.5D0*ULEFT(2)**2/R0)
      P1 = G0*(URIGHT(3)-0.5D0*URIGHT(2)**2/R1)
C     Check for negative pressures ..
      IF (P0.LT.0.0D0 .OR. P1.LT.0.0D0) GO TO 240
      C0 = SQRT(GAMMA*P0/R0)
      C1 = SQRT(GAMMA*P1/R1)
      Z0 = LOG(P0*R0**(-GAMMA))
      Z1 = LOG(P1*R1**(-GAMMA))
C     Calculate PHI1, PHI2 and alpha ..
      IF (PTH.EQ.'P') THEN
         PHI0 = U0 + 2.0D0*G1*C0
         PHI1 = U1 - 2.0D0*G1*C1
      ELSE
         PHI0 = U0 - 2.0D0*G1*C0
         PHI1 = U1 + 2.0D0*G1*C1
      END IF
      AL = EXP(G2*(Z1-Z0))
C     Calculate intermediate values ..
      U13 = (PHI1+AL*PHI0)/(1.0D0+AL)
      U23 = U13
      IF (PTH.EQ.'P') THEN
         C13 = 0.5D0*G0*(PHI0-PHI1)/(1.0D0+AL)
      ELSE
         C13 = -0.5D0*G0*(PHI0-PHI1)/(1.0D0+AL)
      END IF
      C23 = AL*C13
C     Calculate eigenvalues ..
      IF (PTH.EQ.'P') THEN
         D(1) = U0 - C0
         D(2) = U13 - C13
         D(3) = U13
         D(4) = U23 + C23
         D(5) = U1 + C1
      ELSE
         D(1) = U0 + C0
         D(2) = U13 + C13
         D(3) = U13
         D(4) = U23 - C23
         D(5) = U1 - C1
      END IF
      QSP1 = D(1)*D(2)
      QSP2 = D(4)*D(5)
      DO 20 I = 1, 5
         D(I) = SIGN(1.0D0,D(I))
   20 CONTINUE
C
C     Calculate flux at every point ..
      IF (D(1).EQ.-1.0D0) THEN
         DO 40 I = 1, 3
            F(I,1) = 0.0D0
   40    CONTINUE
      ELSE
         W(1,1) = R0
         W(2,1) = U0
         W(3,1) = P0
         F(1,1) = W(1,1)*W(2,1)
         F(2,1) = W(1,1)*W(2,1)**2 + W(3,1)
         F(3,1) = W(2,1)*(GAMMA*W(3,1)*G1+0.5D0*W(1,1)*W(2,1)**2)
      END IF
      IF (D(1).EQ.D(2)) THEN
         DO 60 I = 1, 3
            F(I,2) = 0.0D0
   60    CONTINUE
      ELSE
C       Check for first sonic point ..
         IF (QSP1.LE.0.0D0) THEN
            IF (PTH.EQ.'P') THEN
               CS1 = G3*PHI0
               US1 = CS1
            ELSE
               CS1 = G3*PHI0
               US1 = -CS1
            END IF
            W(1,2) = (CS1**2/(GAMMA*EXP(Z0)))**G1
            W(2,2) = US1
            W(3,2) = EXP(Z0)*W(1,2)**GAMMA
         ELSE
            W(1,2) = 0.0D0
            W(2,2) = 0.0D0
            W(3,2) = 0.0D0
         END IF
         F(1,2) = W(1,2)*W(2,2)
         F(2,2) = W(1,2)*W(2,2)**2 + W(3,2)
         F(3,2) = W(2,2)*(GAMMA*W(3,2)*G1+0.5D0*W(1,2)*W(2,2)**2)
      END IF
      IF (D(2).EQ.D(3)) THEN
         DO 80 I = 1, 3
            F(I,3) = 0.0D0
   80    CONTINUE
      ELSE
         W(1,3) = (C13**2/(GAMMA*EXP(Z0)))**G1
         W(2,3) = U13
         W(3,3) = EXP(Z0)*W(1,3)**GAMMA
         F(1,3) = W(1,3)*W(2,3)
         F(2,3) = W(1,3)*W(2,3)**2 + W(3,3)
         F(3,3) = W(2,3)*(GAMMA*W(3,3)*G1+0.5D0*W(1,3)*W(2,3)**2)
      END IF
      IF (D(3).EQ.D(4)) THEN
         DO 100 I = 1, 3
            F(I,4) = 0.0D0
  100    CONTINUE
      ELSE
         W(1,4) = (C23**2/(GAMMA*EXP(Z1)))**G1
         W(2,4) = U13
         W(3,4) = EXP(Z1)*W(1,4)**GAMMA
         F(1,4) = W(1,4)*W(2,4)
         F(2,4) = W(1,4)*W(2,4)**2 + W(3,4)
         F(3,4) = W(2,4)*(GAMMA*W(3,4)*G1+0.5D0*W(1,4)*W(2,4)**2)
      END IF
      IF (D(4).EQ.D(5)) THEN
         DO 120 I = 1, 3
            F(I,5) = 0.0D0
  120    CONTINUE
      ELSE
C       Check for second sonic point ..
         IF (QSP2.LE.0.0D0) THEN
            IF (PTH.EQ.'P') THEN
               CS2 = G3*PHI1
               US2 = -CS2
            ELSE
               CS2 = G3*PHI1
               US2 = CS2
            END IF
            W(1,5) = (CS2**2/(GAMMA*EXP(Z1)))**G1
            W(2,5) = US2
            W(3,5) = EXP(Z1)*W(1,5)**GAMMA
         ELSE
            W(1,5) = 0.0D0
            W(2,5) = 0.0D0
            W(3,5) = 0.0D0
         END IF
         F(1,5) = W(1,5)*W(2,5)
         F(2,5) = W(1,5)*W(2,5)**2 + W(3,5)
         F(3,5) = W(2,5)*(GAMMA*W(3,5)*G1+0.5D0*W(1,5)*W(2,5)**2)
      END IF
      IF (D(5).EQ.1.0D0) THEN
         DO 140 I = 1, 3
            F(I,6) = 0.0D0
  140    CONTINUE
      ELSE
         W(1,6) = R1
         W(2,6) = U1
         W(3,6) = P1
         F(1,6) = W(1,6)*W(2,6)
         F(2,6) = W(1,6)*W(2,6)**2 + W(3,6)
         F(3,6) = W(2,6)*(GAMMA*W(3,6)*G1+0.5D0*W(1,6)*W(2,6)**2)
      END IF
C     Sum the flux terms ..
      DO 160 I = 1, 3
         IF (PTH.EQ.'P') THEN
            FLUX(I) = F(I,1) + F(I,2)*D(2) + F(I,3) + F(I,4) + F(I,5)
     *                *D(5) + F(I,6)
         ELSE
            FLUX(I) = F(I,1) + F(I,2)*D(2) - F(I,3) - F(I,4) + F(I,5)
     *                *D(5) + F(I,6)
         END IF
  160 CONTINUE
C
      GO TO 260
C
  180 CONTINUE
      ERRMSG =
     *   ' Routine entered with PATH not recognised as ''O'' or ''P''. '
      CALL D02NNQ(ERRMSG,1,0,0,0,0,0.0D0,0.0D0)
      IFAIL1 = 1
      GO TO 260
C
  200 CONTINUE
      ERRMSG =
     *' Routine entered with GAMMA (=R1) less than or equal     to 0.0.
     *'
      CALL D02NNQ(ERRMSG,1,0,0,0,1,GAMMA,0.0D0)
      IFAIL1 = 1
      GO TO 260
C
  220 CONTINUE
      ERRMSG =
     *' Routine entered with left density value i.e. ULEFT(1)   (=R1) an
     *d/or right density value i.e. URIGHT(1) (=R2) less than   0.0. '
      CALL D02NNQ(ERRMSG,1,0,0,0,2,R0,R1)
      IFAIL1 = 2
      GO TO 260
C
  240 CONTINUE
      ERRMSG =
     *' Routine entered with left pressure value (=R1) and/or   right pr
     *essure value (=R2) less than 0.0. '
      CALL D02NNQ(ERRMSG,1,0,0,0,2,P0,P1)
      IFAIL1 = 2
C
  260 CONTINUE
      IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,0,REC)
      RETURN
      END
