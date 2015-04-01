      SUBROUTINE D02RAV(M,N,P,R,ALPHA,A1,B1,X,Y,IY,A2,C2,DEL,FCN,G,
     *                  FCNEP,FCNA,FCNB,JACOBE,JACOBG,A10,B10,GAM,A20,
     *                  B20,JERROR,EPS,IR,IC,UU,RES,MTNMAX,NMAX,MMAX2,F,
     *                  HX,SK,GRADF,AUX,ICA,XAU,LP,MP,LIN,EPSNU,NU,INWT,
     *                  CASI,HMAX,IG,NIG,H,IRN,IP,NIP)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS, EPSNU
      INTEGER           INWT, IY, JERROR, LIN, LP, M, MMAX2, MP, MTNMAX,
     *                  N, NIG, NIP, NMAX, NU, P, R
      LOGICAL           CASI
C     .. Array Arguments ..
      DOUBLE PRECISION  A1(M,M), A10(M,M), A2(MTNMAX,M), A20(M,2),
     *                  ALPHA(M), AUX(MMAX2,M), B1(M,M), B10(M),
     *                  B20(M,2), C2(MTNMAX,M), DEL(M,MTNMAX),
     *                  F(MTNMAX), GAM(M), GRADF(MTNMAX), H(M), HMAX(M),
     *                  HX(NMAX), RES(MTNMAX), SK(MTNMAX), UU(MTNMAX),
     *                  X(NMAX), XAU(MTNMAX), Y(IY,N)
      INTEGER           IC(NMAX,M), ICA(M), IG(NIG), IP(NIP),
     *                  IR(NMAX,M), IRN(M,M)
C     .. Subroutine Arguments ..
      EXTERNAL          FCN, FCNA, FCNB, FCNEP, G, JACOBE, JACOBG
C     .. Local Scalars ..
      DOUBLE PRECISION  DT, GNOR, PNOR, RABS, RABS1, REOLD, RMAX, ROL,
     *                  ROL1, S, SCPR, T, TMIN, TN
      INTEGER           I, I1, J, J1, K1, K1J, K1JM, K1JMMP, MPN, MPNM,
     *                  MPNMPJ, NADV, P1
      LOGICAL           SING
C     .. Local Arrays ..
      CHARACTER*80      REC(5)
C     .. External Subroutines ..
      EXTERNAL          D02RAU, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Executable Statements ..
C
C     ERROR EXIT
C     JERROR = 3    NON  CONVERGENCE
C
C     ***        NEWTON ITERATION STARTS        ***
C
      DT = 1.D0
      P1 = P + 1
      T = 1.D0
      TMIN = 2.D0**(-8)
      SCPR = 0.D0
      MPN = M*N
      DO 20 I = 1, MPN
         UU(I) = 0.D0
   20 CONTINUE
      MPNM = M*(N-1)
      INWT = 0
      ROL = 1.D20
      REOLD = 1.0D20
      IF (MP.EQ.0) GO TO 40
      CALL X04ABF(0,NADV)
      IF (LIN.NE.1) THEN
         CALL X04BAF(NADV,' ')
         WRITE (REC,FMT=99999)
         CALL X04BAF(NADV,REC(1))
      END IF
      IF (LIN.GE.3) THEN
         WRITE (REC,FMT=99998) N
         CALL X04BAF(NADV,REC(1))
      END IF
      IF (LIN.NE.1) THEN
         WRITE (REC,FMT=99997) NU, EPS
         CALL X04BAF(NADV,REC(1))
      END IF
C
C     RESIDUAL  COMPUTATION
C
   40 RABS = 0.D0
      GO TO (180,100,60,60) LIN
   60 DO 80 I = 1, N
         I1 = (I-1)*M + 1
         CALL FCNEP(X(I),EPSNU,Y(1,I),F(I1),M)
   80 CONTINUE
      CALL G(EPSNU,Y(1,1),Y(1,N),ALPHA,M)
      GO TO 300
  100 DO 120 I = 1, N
         I1 = (I-1)*M + 1
         CALL FCN(X(I),Y(1,I),F(I1))
  120 CONTINUE
      J = 0
      DO 140 I = 1, M
         IF (B20(I,1).NE.0.0D0) GO TO 140
         J = J + 1
         ALPHA(J) = Y(I,1) - A20(I,1)
  140 CONTINUE
      DO 160 I = 1, M
         IF (B20(I,2).NE.0.0D0) GO TO 160
         J = J + 1
         ALPHA(J) = Y(I,N) - A20(I,2)
  160 CONTINUE
      GO TO 300
  180 DO 240 I = 1, N
         CALL FCNA(X(I),A10)
         CALL FCNB(X(I),B10)
         DO 220 J = 1, M
            S = 0.0D0
            DO 200 J1 = 1, M
               S = S + A10(J,J1)*Y(J1,I)
  200       CONTINUE
            I1 = (I-1)*M + J
            F(I1) = S + B10(J)
  220    CONTINUE
  240 CONTINUE
      DO 280 I = 1, M
         S = 0.0D0
         DO 260 J = 1, M
            S = S + A1(I,J)*Y(J,1) + B1(I,J)*Y(J,N)
  260    CONTINUE
         ALPHA(I) = S - GAM(I)
  280 CONTINUE
  300 CONTINUE
      RMAX = 0.0D0
      IF (P.EQ.0) GO TO 360
      DO 320 I = 1, P
         RES(I) = -ALPHA(I)
         RMAX = MAX(RMAX,ABS(RES(I)))
  320 CONTINUE
      IF (RMAX.EQ.0.0D0) GO TO 360
      RABS1 = 0.0D0
      DO 340 I = 1, P
         RABS1 = RABS1 + (RES(I)/RMAX)**2
  340 CONTINUE
      RABS = RABS + (RABS1*RMAX)*RMAX
C
  360 DO 420 I = 2, N
         I1 = I - 1
         K1 = I1*M
         RMAX = 0.0D0
         DO 380 J = 1, M
            K1J = K1 + J
            K1JM = K1J - M
            K1JMMP = K1J - M + P
            RES(K1JMMP) = -Y(J,I) + Y(J,I-1) + .5D0*HX(I1)*(F(K1J)
     *                    +F(K1JM)) + SK(K1JM)
            RMAX = MAX(RMAX,ABS(RES(K1JMMP)))
  380    CONTINUE
         IF (RMAX.EQ.0.0D0) GO TO 420
         RABS1 = 0.0D0
         DO 400 J = 1, M
            K1JMMP = K1 + J - M + P
            RABS1 = RABS1 + (RES(K1JMMP)/RMAX)**2
  400    CONTINUE
         RABS = RABS + (RABS1*RMAX)*RMAX
  420 CONTINUE
C
      RMAX = 0.0D0
      DO 440 J = P1, M
         MPNMPJ = MPNM + J
         RES(MPNMPJ) = -ALPHA(J)
         RMAX = MAX(RMAX,ABS(ALPHA(J)))
  440 CONTINUE
      IF (RMAX.EQ.0.0D0) GO TO 480
      RABS1 = 0.0D0
      DO 460 J = P1, M
         RABS1 = RABS1 + (ALPHA(J)/RMAX)**2
  460 CONTINUE
      RABS = RABS + (RABS1*RMAX)*RMAX
  480 RABS1 = SQRT(RABS)
      IF (MP.NE.0 .AND. LIN.NE.1) THEN
         WRITE (REC,FMT=99996) INWT, RABS1
         CALL X04BAF(NADV,REC(1))
      END IF
      IF (INWT.EQ.0 .AND. (EPSNU.EQ.1.D0 .OR. EPSNU.EQ.0.D0))
     *    GO TO 700
C
C     CHECK FOR CONVERGENCE
C
      IF (RABS1.GT.EPS) GO TO 500
C     JACOBIAN IS KEPT CONSTANT UNTIL NEXT MESH CHANGE.
C     STEP AND ANGLE CONTROL ARE SHORTCIRCUITED.
      IF (EPSNU.GE.1.D0) CASI = .TRUE.
      GO TO 1020
C
C
C     .....                       NEWTON EXIT
C     CONVERGENCE  OR  TOO  MANY  ITERATIONS ...
C
C
C     NEWTON TEST IN ORDER TO AVOID CYCLING
C
  500 IF ((REOLD-RABS.GE..5D0*T*SCPR .AND. INWT.LT.15)
     *    .OR. (NU.GT.0 .AND. INWT.EQ.1)) GO TO 700
      IF (INWT.EQ.15) GO TO 580
C
C     STEP CONTROL STARTS ****
C
      IF ( .NOT. (CASI .OR. LIN.EQ.1)) GO TO 520
      IF (RABS1.LE.100.D0*EPS) GO TO 600
      GO TO 680
  520 TN = .5D0*T
      IF (TN.LT.TMIN .OR. RABS.GT.ROL) GO TO 580
      ROL = RABS
      DT = TN - T
      T = TN
      I1 = 0
      DO 560 I = 1, N
         DO 540 J = 1, M
            I1 = I1 + 1
            Y(J,I) = Y(J,I) + DT*UU(I1)
  540    CONTINUE
  560 CONTINUE
C     IF (MP.NE.0 .AND. LIN.GE.3) THEN
C        WRITE (REC,99995) T, DT
C        CALL X04BAF(NADV,REC(1))
C     ENDIF
      GO TO 40
C
  580 ROL1 = SQRT(ROL)
      IF (ROL1.GE.100.D0*EPS) GO TO 680
      RABS1 = ROL1
  600 I1 = 0
      DO 640 I = 1, N
         DO 620 J = 1, M
            I1 = I1 + 1
            Y(J,I) = Y(J,I) - DT*UU(I1)
  620    CONTINUE
  640 CONTINUE
      IF (MP.NE.0) THEN
         CALL X04BAF(NADV,' ')
         WRITE (REC,FMT=99995) NU, RABS1, EPS
         DO 660 I = 1, 5
            CALL X04BAF(NADV,REC(I))
  660    CONTINUE
      END IF
C     NEWTON DID NOT QUITE REACH THE TOLERANCE.
C     FURTHER MESH REFINEMENTS ARE NOT ALLOWED.
      JERROR = 4
C     JACOBIAN IS KEPT CONSTANT UNTIL NEXT MESH CHANGE.
C     STEP AND ANGLE CONTROL ARE SHORTCIRCUITED.
      IF (EPSNU.GE.1.D0) CASI = .TRUE.
      GO TO 1020
C     *****
C
C     WE ASSUME DIVERGENCE AND RETURN
C
C     ...    ERROR  EXIT  3    ........
C
  680 JERROR = 3
      GO TO 1020
C
  700 CALL D02RAU(M,N,P,R,X,Y,IY,FCN,G,FCNEP,JACOBG,JACOBE,B20,FCNA,A1,
     *            B1,A2,C2,DEL,CASI,SING,IR,IC,UU,RES,LIN,MTNMAX,NMAX,
     *            MMAX2,HX,GRADF,AUX,ICA,XAU,EPSNU,HMAX,IG,NIG,H,IRN,IP,
     *            NIP,F,ALPHA,JERROR,LP)
      IF (JERROR.EQ.6 .OR. JERROR.EQ.8) RETURN
      IF (SING .AND. MP.NE.0) THEN
         WRITE (REC,FMT=99994)
         CALL X04BAF(NADV,REC(1))
      END IF
      REOLD = RABS
      SCPR = 1.D-10
      T = 1.D0
      IF (CASI .OR. LIN.EQ.1) GO TO 920
      GNOR = 0.D0
      SCPR = 0.D0
      PNOR = 0.D0
      RMAX = 0.0D0
      DO 720 I = 1, MPN
         RMAX = MAX(RMAX,ABS(UU(I)))
  720 CONTINUE
      IF (RMAX.EQ.0.0D0) GO TO 760
      DO 740 I = 1, MPN
         PNOR = PNOR + (UU(I)/RMAX)**2
  740 CONTINUE
      PNOR = (PNOR*RMAX)*RMAX
  760 RMAX = 0.0D0
      DO 780 I = 1, MPN
         RMAX = MAX(RMAX,ABS(GRADF(I)))
  780 CONTINUE
      IF (RMAX.EQ.0.0D0) GO TO 820
      DO 800 I = 1, MPN
         GNOR = GNOR + (GRADF(I)/RMAX)**2
  800 CONTINUE
      GNOR = (GNOR*RMAX)*RMAX
  820 DO 840 I = 1, MPN
         SCPR = SCPR + GRADF(I)*UU(I)
  840 CONTINUE
      IF (MP.NE.0 .AND. LIN.GE.3) THEN
         WRITE (REC,FMT=99993) PNOR, GNOR, SCPR
         DO 860 I = 1, 3
            CALL X04BAF(NADV,REC(I))
  860    CONTINUE
      END IF
      IF (PNOR.EQ.0.D0) GO TO 580
C
C     WE CHECK IF THE DIRECTION  UU  IS OF DESCENT (AS IT SHOULD)
C     AND ALSO IF THE IDENTITY   GRADF,UU =RABS  IS APPROXIMATELY
C     VERIFIED. IF EITHER ONE OF THESE CHECKS FAIL, WE TAKE  GRADF
C     AS OUR NEXT SEARCH DIRECTION, SINCE THE ABOVE INDICATES THAT
C     UU  IS AN UNRELIABLE DIRECTION DUE TO ILL-CONDITIONING IN
C     THE JACOBIAN.
C
      IF (SCPR) 900, 900, 880
  880 IF (ABS(SCPR-RABS).LE..1D0*RABS) GO TO 920
  900 SCPR = -1.D0
C
C     APPROXIMATE SOLUTION IS CORRECTED
C
  920 IF (SCPR.GT.0.D0) GO TO 960
      DO 940 I = 1, MPN
         UU(I) = GRADF(I)
  940 CONTINUE
  960 I1 = 0
      DO 1000 I = 1, N
         DO 980 J = 1, M
            I1 = I1 + 1
            Y(J,I) = Y(J,I) + UU(I1)
  980    CONTINUE
 1000 CONTINUE
      IF (SCPR.LE.0.D0) SCPR = GNOR
      INWT = INWT + 1
      GO TO 40
 1020 RETURN
C
C     9995 FORMAT (23H     STEP CONTROL   T =, 1P,E10.2, 7H   DT =,
C     * 1PE10.2)
99999 FORMAT ('  MONITORING NEWTON ITERATION')
99998 FORMAT ('   NUMBER OF POINTS IN CURRENT MESH =',I5)
99997 FORMAT ('   CORRECTION NUMBER ',I4,'   RESIDUAL SHOULD BE .LE.',
     *  1P,D10.2)
99996 FORMAT ('    ITERATION NUMBER ',I4,'   RESIDUAL =',1P,D10.2)
99995 FORMAT ('   WARNING --- CORRECTION NUMBER ',I4,/'    RESIDUAL CO',
     *  'ULD ONLY BE REDUCED TO ',1P,D10.2,/'    INSTEAD OF BELOW ',1P,
     *  D10.2,/'    IF ON TERMINATION TOL HAS NOT BEEN REACHED',/'    ',
     *  'THE PROBLEM REQUIRES A HIGHER MACHINE PRECISION')
99994 FORMAT ('   SINGULAR JACOBIAN')
99993 FORMAT ('     SQUARED NORM OF CORRECTION =',1P,D10.2,/'     SQUA',
     *  'RED NORM OF GRADIENT   =',1P,D10.2,/'     SCALAR PRODUCT OF C',
     *  'ORRECTION AND GRADIENT =',1P,D10.2)
      END
