      SUBROUTINE D02RAY(M,NMAX,N,P,R,MTNMAX,MMAX2,A,B,TOL,X,Y,IY,ABT,
     *                  FCN,G,FCNEP,FCNA,FCNB,JACOBE,JACOBG,JACEPS,
     *                  JACGEP,A10,B10,GAM,A20,B20,ALPHA,A1,B1,EJ,A2,C2,
     *                  DEL,UU,RES,F,HX,SK,GRADF,AUX,XAU,IC,IR,IQJ,ICA,
     *                  DELEPS,LP,MP,INIT,LIN,H2,HMAX,IG,NIG,IP2,NIP,
     *                  IRN,IFLAG)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 9B REVISED. IER-360 (JAN 1982)
C     MARK 9C REVISED. IER-368 (MAY 1982).
C     MARK 9C REVISED. IER-373 (JUN 1982).
C     MARK 10 REVISED. IER-377 (JUN 1982).
C     MARK 10 REVISED. IER-378 (JUN 1982).
C     MARK 10B REVISED. IER-401 (JAN 1983).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     FCN, FCNA, FCNB, FCNEP, G, JACEPS, JACGEP, JACOBE, JACOBG
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, DELEPS, TOL
      INTEGER           IFLAG, INIT, IY, LIN, LP, M, MMAX2, MP, MTNMAX,
     *                  N, NIG, NIP, NMAX, P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A1(M,M), A10(M,M), A2(MTNMAX,M), A20(M,2),
     *                  ABT(M), ALPHA(M), AUX(MMAX2,M), B1(M,M), B10(M),
     *                  B20(M,2), C2(MTNMAX,M), DEL(M,MTNMAX), EJ(NMAX),
     *                  F(MTNMAX), GAM(M), GRADF(MTNMAX), H2(M),
     *                  HMAX(M), HX(NMAX), RES(MTNMAX), SK(MTNMAX),
     *                  UU(MTNMAX), X(NMAX), XAU(MTNMAX), Y(IY,NMAX)
      INTEGER           IC(NMAX,M), ICA(M), IG(NIG), IP2(NIP),
     *                  IQJ(NMAX), IR(NMAX,M), IRN(M,M)
C     .. Subroutine Arguments ..
      EXTERNAL          FCN, FCNA, FCNB, FCNEP, G, JACEPS, JACGEP,
     *                  JACOBE, JACOBG
C     .. Local Scalars ..
      DOUBLE PRECISION  ALG, AUXI, BMA, C3, E, EPBAR, EPS, EPS1, EPSMAC,
     *                  EPSMCH, EPSNU, ERRNEW, ERROLD, H, HI, ONE, RABS,
     *                  S, SIG02, TE, TEM, U1, UUN, XN, XTE, XXN, YKI,
     *                  ZERO
      INTEGER           I, I1, ICON, IERROR, IFAIL, IFI, IFIN, II, IIM,
     *                  IIMPJ, IIP, IIPJ, IIPPJ, IMNETC, IND, INWT, IP,
     *                  IPM, IPP, IQ, ITEMP, J, J1, JI, JJ, K, K11, KI,
     *                  KII1, KII1PK, KIJ, KK, KKPJ, KMAX, L, MPN, MPNM,
     *                  MPNMPI, N1, N1TMPI, N2, NADV, NERR, NMA, NOLD,
     *                  NP, NPU, NTEM, NTOP, NU, NU1, NU2
      LOGICAL           CASI, SING
C     .. Local Arrays ..
      DOUBLE PRECISION  AA(50,5), EPSNU1(1), H1(1), HMAX1(1)
      INTEGER           IG1(3), IP1(2)
      CHARACTER*80      REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          D02RAR, D02RAS, D02RAU, D02RAV, D02RAW, D02RAX,
     *                  X04AAF, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, DBLE, SQRT
C     .. Data statements ..
      DATA              ZERO/0.0D0/, ONE/1.0D0/
C     .. Executable Statements ..
      IF (LP.NE.0) CALL X04AAF(0,NERR)
      IF (MP.EQ.0) GO TO 20
      CALL X04ABF(0,NADV)
      CALL X04BAF(NADV,' ')
      IF (LIN.GE.3) WRITE (REC,FMT=99985)
      IF (LIN.EQ.1) WRITE (REC,FMT=99984)
      IF (LIN.EQ.2) WRITE (REC,FMT=99983)
      CALL X04BAF(NADV,REC(1))
   20 EPSMCH = X02AJF()
C
C     INITIALIZATION
C
      CASI = .FALSE.
      SING = .FALSE.
      MPN = M*N
C
C     ***********************************************************
C     THIS CONSTANT IS MACHINE DEPENDENT=
C     EPSMAC  IS APPROXIMATELY 10* RELATIVE PRECISION OF
C     FLOATING POINT ARITHMETIC
      EPSMAC = 10.0D0*EPSMCH
C     ***********************************************************
C
      EPSMAC = MAX(EPSMAC,1.D-4*TOL)
      N1 = N - 1
      EPSNU = ZERO
      IF (DELEPS.EQ.ZERO) EPSNU = ONE
      RABS = ONE
      IFLAG = 0
      EPBAR = 1.D10
      ICON = 3
      NOLD = 1
      NMA = NMAX
      NTOP = (NMAX*9)/10
      NTOP = MAX(NTOP,NMAX-16)
      BMA = 0.07D0
C     INITIALIZATION OF ARRAYS USED IN DEFFERRED CORRECTION.
      AA(1,2) = ONE
      DO 40 I = 2, 50
         AA(I,2) = ZERO
   40 CONTINUE
      DO 60 NU1 = 1, 20
         NU2 = 2*NU1 + 1
         AA(NU2,1) = -(DBLE(NU1)/DBLE(NU2))*0.5D0**(2*NU1-1)
         AA(NU2+1,1) = ZERO
   60 CONTINUE
      IF (INIT.NE.0 .AND. LIN.EQ.1) GO TO 100
      IF (INIT.NE.0) GO TO 160
C
C     FIRST  APPROXIMATION  FOR  Y  AND  X (INIT .EQ. 0)
C
      X(1) = A
      X(N) = B
      H = (B-A)/DBLE(N1)
      DO 80 I = 2, N1
         X(I) = A + DBLE(I-1)*H
   80 CONTINUE
  100 DO 140 I = 1, N
         DO 120 J = 1, M
            Y(J,I) = ZERO
  120    CONTINUE
  140 CONTINUE
C
  160 H = ZERO
      DO 180 I = 1, N1
         HX(I) = X(I+1) - X(I)
         IF (HX(I).GT.H) H = HX(I)
  180 CONTINUE
C     INITIALISATION OF ARRAYS IN JACOBIAN CALCUATION
      IF (LIN.EQ.1 .OR. LIN.EQ.3) GO TO 240
      DO 220 I = 1, M
         IP2(I) = (I-1)*M + 1
         HMAX(I) = MAX(ABS(Y(I,1)),1.0D0)
         H2(I) = SQRT(EPSMCH)*HMAX(I)
         DO 200 J = 1, M
            IRN(I,J) = I
  200    CONTINUE
  220 CONTINUE
      IG(1) = 0
      IP2(M+1) = M*M + 1
      IF (DELEPS.EQ.0.0D0 .OR. LIN.EQ.2) GO TO 240
      IG1(1) = 0
      IP1(1) = 1
      IP1(2) = M + 1
      HMAX1(1) = 1.0D0
      H1(1) = SQRT(EPSMCH)
C
C     ....  MAIN BODY OF D02RAY ......
C
C     NU IS CURRENT NUMBER OF CORRECTIONS.
  240 NU = 0
      EPS = MAX(EPSMAC,H**2*.1D0)
C
C     ....      ENTER  AFTER  MESH  CHANGE    .......
C
  260 ERROLD = 1.0D20
C     KMAX IS MAXIMUM NUMBER OF CORRECTIONS.
      KMAX = 20
      MPN = M*N
      N1 = N - 1
      MPNM = M*N1
      IPM = MPNM + P + 1
      C3 = 0.8D0
      DO 280 I = 1, MPN
         SK(I) = ZERO
  280 CONTINUE
      IF (NU.EQ.0) GO TO 600
C
C     AFTER  MESH  CHANGE  WE HAVE TO INITIALIZE  SK
C     IF  NU .GT. 0
C
      GO TO (420,340,300,300) LIN
  300 DO 320 I = 1, N
         I1 = (I-1)*M + 1
         CALL FCNEP(X(I),EPSNU,Y(1,I),F(I1),M)
  320 CONTINUE
      CALL G(EPSNU,Y(1,1),Y(1,N),ALPHA,M)
      GO TO 540
  340 DO 360 I = 1, N
         I1 = (I-1)*M + 1
         CALL FCN(X(I),Y(1,I),F(I1))
  360 CONTINUE
      J = 0
      DO 380 I = 1, M
         IF (B20(I,1).NE.0.0D0) GO TO 380
         J = J + 1
         ALPHA(J) = Y(I,1) - A20(I,1)
  380 CONTINUE
      DO 400 I = 1, M
         IF (B20(I,2).NE.0.0D0) GO TO 400
         J = J + 1
         ALPHA(J) = Y(I,N) - A20(I,2)
  400 CONTINUE
      GO TO 540
  420 DO 480 I = 1, N
         CALL FCNA(X(I),A10)
         CALL FCNB(X(I),B10)
         DO 460 J = 1, M
            S = 0.0D0
            DO 440 J1 = 1, M
               S = S + A10(J,J1)*Y(J1,I)
  440       CONTINUE
            I1 = (I-1)*M + J
            F(I1) = S + B10(J)
  460    CONTINUE
  480 CONTINUE
      DO 520 I = 1, M
         S = 0.0D0
         DO 500 J = 1, M
            S = S + A1(I,J)*Y(J,1) + B1(I,J)*Y(J,N)
  500    CONTINUE
         ALPHA(I) = S - GAM(I)
  520 CONTINUE
  540 CONTINUE
      CALL D02RAX(NU,2,2,N1,M,AA,X,NMAX,F,RES,MTNMAX,AA(1,4),AA(1,5)
     *            ,IERROR)
      DO 580 I = 1, N1
         KI = (I-1)*M
         DO 560 J = 1, M
            KIJ = KI + J
            SK(KIJ) = HX(I)*RES(KIJ)
  560    CONTINUE
  580 CONTINUE
      IF (NU.LT.KMAX) GO TO 600
      NU = NU - 1
      GO TO 1420
C
C     NU IS TOO LARGE  GO TO REFINE THE MESH ...
C
C     ***  NEWTON  ITERATION  ***
C
  600 IF (EPSNU.GE.ONE) GO TO 640
      EPS1 = EPS
  620 EPS = MAX(RABS,EPS)
  640 CALL D02RAV(M,N,P,R,ALPHA,A1,B1,X,Y,IY,A2,C2,DEL,FCN,G,FCNEP,FCNA,
     *            FCNB,JACOBE,JACOBG,A10,B10,GAM,A20,B20,IFLAG,EPS,IR,
     *            IC,UU,RES,MTNMAX,NMAX,MMAX2,F,HX,SK,GRADF,AUX,ICA,XAU,
     *            LP,MP,LIN,EPSNU,NU,INWT,CASI,HMAX,IG,NIG,H2,IRN,IP2,
     *            NIP)
      IF (IFLAG.NE.3 .AND. IFLAG.NE.6 .AND. IFLAG.NE.8) GO TO 700
C     NEWTON STEP FAILS TO CONVERGE.
      IF (LP.EQ.0 .OR. IFLAG.NE.3) RETURN
      CALL X04BAF(NERR,' ')
      WRITE (REC,FMT=99999)
      DO 660 I = 1, 2
         CALL X04BAF(NERR,REC(I))
  660 CONTINUE
      IF (DELEPS.EQ.0.0D0) THEN
         WRITE (REC,FMT=99998)
         CALL X04BAF(NERR,REC(1))
      END IF
      IF (MP.EQ.0) THEN
         WRITE (REC,FMT=99997)
         DO 680 I = 1, 2
            CALL X04BAF(NERR,REC(I))
  680    CONTINUE
      END IF
      RETURN
  700 IF (EPSNU.GE.ONE) GO TO 1120
C
C     CONTINUATION  -  CHOOSE STEP AND NEW INITIAL PROFILE
C
      IF (LIN.EQ.4) GO TO 740
      DO 720 I = 1, N
         I1 = (I-1)*M + 1
         CALL JACEPS(X(I),EPSNU,Y(1,I),F(I1),M)
  720 CONTINUE
      CALL JACGEP(EPSNU,Y(1,1),Y(1,N),ALPHA,M)
      GO TO 860
  740 EPSNU1(1) = EPSNU
      DO 800 I = 1, N
         IFAIL = 1
         I1 = (I-1)*M + 1
         CALL FCNEP(X(I),EPSNU1(1),Y(1,I),ABT,M)
         IND = 1
  760    CALL D02RAR(M,1,IRN,M,IP1,2,H1,EPSNU1,ABT,HMAX1,10.0D0,100.0D0,
     *               1000.0D0,EPSMCH,EPSMCH,.TRUE.,F(I1)
     *               ,IG1,3,AUX,M+1,XAU,IND,IFAIL)
         IF (IND.EQ.0) GO TO 780
         CALL FCNEP(X(I),EPSNU1(1),Y(1,I),AUX,M)
         GO TO 760
  780    IF (IFAIL.EQ.0) GO TO 800
         IFLAG = 6
         RETURN
  800 CONTINUE
      IND = 1
      I1 = (N-1)*M + 1
      CALL G(EPSNU1(1),Y(1,1),Y(1,N),ABT,M)
      IFAIL = 1
  820 CALL D02RAR(M,1,IRN,M,IP1,2,H1,EPSNU1,ABT,HMAX1,10.0D0,100.0D0,
     *            1000.0D0,EPSMCH,EPSMCH,.TRUE.,ALPHA,IG1,3,AUX,M+1,XAU,
     *            IND,IFAIL)
      IF (IND.EQ.0) GO TO 840
      CALL G(EPSNU1(1),Y(1,1),Y(1,N),AUX,M)
      GO TO 820
  840 IF (IFAIL.EQ.0) GO TO 860
      IFLAG = 6
      RETURN
  860 CONTINUE
      DO 880 I = 1, P
         UU(I) = -ALPHA(I)
  880 CONTINUE
      DO 920 I = 1, N1
         II = (I-1)*M
         IIP = II + P
         IIM = II + M
         DO 900 J = 1, M
            IIPPJ = IIP + J
            IIMPJ = IIM + J
            IIPJ = II + J
            UU(IIPPJ) = .5D0*HX(I)*(F(IIMPJ)+F(IIPJ))
  900    CONTINUE
  920 CONTINUE
      IP = P + 1
      DO 940 I = IP, M
         N1TMPI = N1*M + I
         UU(N1TMPI) = -ALPHA(I)
  940 CONTINUE
      J = N*M
      DO 960 I = 1, J
         IF (UU(I).NE.0.0D0) GO TO 1000
  960 CONTINUE
      IF (LP.NE.0) THEN
         CALL X04BAF(NERR,' ')
         WRITE (REC,FMT=99996) EPSNU, DELEPS
         DO 980 I = 1, 2
            CALL X04BAF(NERR,REC(I))
  980    CONTINUE
      END IF
      IFLAG = 7
      RETURN
C     SOLVE  A *UU = UU.
 1000 CALL D02RAS(A2,C2,DEL,UU,M,N,P,R,IR,IC,UU,MTNMAX,NMAX,XAU,SING)
      IF (SING .AND. MP.NE.0) THEN
         WRITE (REC,FMT=99981)
         CALL X04BAF(NADV,REC(1))
      END IF
      IF (INWT-1) 1020, 1020, 1040
 1020 DELEPS = 2.D0*DELEPS
      GO TO 1060
 1040 DELEPS = DELEPS/DBLE(INWT-1)
      IF (DELEPS.GE.SQRT(EPSMCH)) GO TO 1060
      IFLAG = 9
      IF (LP.NE.0) THEN
         CALL X04BAF(NERR,' ')
         WRITE (REC,FMT=99994) DELEPS
         CALL X04BAF(NERR,REC(1))
      END IF
      RETURN
 1060 I1 = 0
      DO 1100 I = 1, N
         DO 1080 J = 1, M
            I1 = I1 + 1
            Y(J,I) = Y(J,I) + DELEPS*UU(I1)
 1080    CONTINUE
 1100 CONTINUE
      EPSNU = MIN(EPSNU+DELEPS,ONE)
      IF (MP.NE.0) THEN
         CALL X04BAF(NADV,' ')
         WRITE (REC,FMT=99995) EPSNU, DELEPS
         CALL X04BAF(NADV,REC(1))
      END IF
      IF (EPSNU.LT.ONE) GO TO 620
      EPS = EPS1
      GO TO 640
C
C     CORRECTION AND ERROR CONTROL STARTS
C
 1120 CALL D02RAX(NU+1,2,2,N1,M,AA,X,NMAX,F,RES,MTNMAX,AA(1,4),AA(1,5)
     *            ,IERROR)
C     JUMP IF WE HAVE EXCEEDED THE NUMBER OF CORRECTIONS ALLOWED.
      IF (IERROR.EQ.1) GO TO 1420
      DO 1160 I = 1, N1
         II = (I-1)*M
         DO 1140 J = 1, M
            KI = II + J
            AUXI = RES(KI)*HX(I)
            RES(KI) = SK(KI) - AUXI
            SK(KI) = AUXI
 1140    CONTINUE
 1160 CONTINUE
      IF (ICON.LE.12) GO TO 1600
 1180 IF (P.EQ.0) GO TO 1220
      DO 1200 I = 1, P
         UU(I) = ZERO
 1200 CONTINUE
 1220 DO 1240 I = IPM, MPN
         UU(I) = ZERO
 1240 CONTINUE
      DO 1260 I = 1, MPNM
         IPP = I + P
         UU(IPP) = RES(I)
 1260 CONTINUE
C     SOLVE   A *UU = UU.
      CALL D02RAU(M,N,P,R,X,Y,IY,FCN,G,FCNEP,JACOBG,JACOBE,B20,FCNA,A1,
     *            B1,A2,C2,DEL,.TRUE.,SING,IR,IC,UU,UU,LIN,MTNMAX,NMAX,
     *            MMAX2,HX,GRADF,AUX,ICA,XAU,EPSNU,HMAX,IG,NIG,H2,IRN,
     *            IP2,NIP,F,ALPHA,IFLAG,LP)
      IF (SING .AND. MP.NE.0) THEN
         WRITE (REC,FMT=99981)
         CALL X04BAF(NADV,REC(1))
      END IF
      IF (IFLAG.EQ.6 .OR. IFLAG.EQ.8) RETURN
C
C     ESTIMATE FOR MAX. ABSOLUTE ERROR (BY COMPONENTS)
C
      ICON = 15
      ERRNEW = ZERO
      DO 1280 J = 1, M
         ABT(J) = ZERO
 1280 CONTINUE
      DO 1320 I = 1, N
         KK = (I-1)*M
         DO 1300 J = 1, M
            KKPJ = KK + J
            U1 = ABS(UU(KKPJ))
            IF (U1.GT.ABT(J)) ABT(J) = U1
 1300    CONTINUE
 1320 CONTINUE
      DO 1340 J = 1, M
         IF (ABT(J).GT.ERRNEW) ERRNEW = ABT(J)
 1340 CONTINUE
      K = NU + 1
      IF (MP.EQ.0) GO TO 1380
      IF (LIN.LE.2) THEN
         WRITE (REC,FMT=99990) N
         CALL X04BAF(NADV,REC(1))
      END IF
      CALL X04BAF(NADV,' ')
      WRITE (REC,FMT=99993) NU, ERRNEW
      CALL X04BAF(NADV,REC(1))
      WRITE (REC,FMT=99992)
      CALL X04BAF(NADV,REC(1))
      DO 1360 I = 1, M, 7
         WRITE (REC,FMT=99991) (ABT(J),J=I,MIN(M,I+6))
         CALL X04BAF(NADV,REC(1))
 1360 CONTINUE
 1380 IF (ERRNEW.LE.TOL) RETURN
C
C     ....           PRECISION  ACHIEVED      ............
C
C     IF NOT ENOUGH POINTS  REFINE THE MESH ******
      IF (NU+1.GE.KMAX) GO TO 1600
      IF (ERRNEW.LE..1D0*ERROLD) GO TO 1400
      IF (ERRNEW.GT.C3*ERROLD) GO TO 1420
      C3 = 0.5D0*C3
C
C     EITHER KEEP CORRECTING ...
C
 1400 ERROLD = ERRNEW
      EPS = MAX(EPSMAC,1.D-3*ERROLD)
      NU = NU + 1
      GO TO 640
C
C     OR REFINE THE MESH, UNLESS IFLAG = 4 ...
C
 1420 IF (IFLAG.NE.4) GO TO 1460
      IF (LP.NE.0) THEN
         CALL X04BAF(NERR,' ')
         WRITE (REC,FMT=99989)
         DO 1440 I = 1, 3
            CALL X04BAF(NERR,REC(I))
 1440    CONTINUE
      END IF
      RETURN
 1460 IF (N.LE.NTOP) GO TO 1480
C
C     ....           TOO MANY GRID POINTS        ......
C
      IFLAG = 2
      IF (LP.NE.0) THEN
         CALL X04BAF(NERR,' ')
         WRITE (REC,FMT=99988)
         CALL X04BAF(NERR,REC(1))
      END IF
      RETURN
 1480 EPS = MAX(EPSMAC,1.D-3*ERROLD)
      IF (ERROLD.LE.ONE) GO TO 1500
      EPBAR = MIN(ONE,1.D-2*ERROLD)
      GO TO 1520
 1500 EPBAR = .01D0*ERROLD
 1520 ICON = 0
      NOLD = N
      BMA = ONE
      IF (NU.LT.1) GO TO 1600
      DO 1540 I = 1, MPN
         SK(I) = SK(I) + RES(I)
 1540 CONTINUE
      NU = NU - 1
      IF (NU.EQ.0) GO TO 1600
      CALL D02RAX(NU,2,2,N1,M,AA,X,NMAX,F,RES,MTNMAX,AA(1,4),AA(1,5)
     *            ,IERROR)
      DO 1580 I = 1, N1
         II = (I-1)*M
         DO 1560 J = 1, M
            KI = II + J
            RES(KI) = RES(KI)*HX(I) - SK(KI)
 1560    CONTINUE
 1580 CONTINUE
C
C     ***** MESH VARIATION *****
C     EQUIDISTRIBUTION OF THE  L2  NORM OF THE ERROR
C     FOR THE O(H**(2*NU+2))  METHOD.
C
 1600 ICON = ICON + 1
      IF (MP.NE.0 .AND. LIN.GE.3) THEN
         CALL X04BAF(NADV,' ')
         WRITE (REC,FMT=99987)
         CALL X04BAF(NADV,REC(1))
      END IF
      ALG = 1.5D0
      SIG02 = ONE/DBLE(2*NU+2)
      TEM = ZERO
      UUN = ZERO
      DO 1640 I = 1, N1
         TE = ZERO
         KI = (I-1)*M
         DO 1620 J = 1, M
            K11 = KI + J
            TEM = MAX(TEM,ABS(Y(J,I)))
            TE = MAX(TE,ABS(RES(K11)))
 1620    CONTINUE
         EJ(I) = (TE/HX(I))**SIG02
         UUN = UUN + EJ(I)
 1640 CONTINUE
C
      IF (ICON.GT.1 .AND. NOLD.GT.1) GO TO 1660
      EPBAR = MAX(MIN(EPBAR,TEM*BMA),TOL)
      E = EPBAR**SIG02
C     IF (MP.NE.0 .AND. LIN.GE.3) WRITE (NADV,99981) E, UUN, TEM,
C     * EPBAR
 1660 IQ = 0
      N2 = N - 2
      I = 0
 1680 I = I + 1
      IQJ(I) = EJ(I)/E - 0.33D0
      IQ = IQ + IQJ(I)
      IF (I.LE.N2) GO TO 1680
      IFIN = (4*N)/100
      NMA = MIN(NMAX-N,70)
      IF (MP.NE.0 .AND. LIN.GE.3) THEN
         WRITE (REC,FMT=99982) IQ
         CALL X04BAF(NADV,REC(1))
      END IF
      IF (IQ.LE.IFIN .OR. NMA.LE.0) GO TO 1920
      IF (IQ.LE.NMA) GO TO 1700
C
C     WE ATTEMPT TO DIMINISH THE NUMBER OF POINTS TO BE
C     INTRODUCED
C
      IF (ALG.LT.1.09D0) GO TO 1920
C     REDUCE ALG AND TRY AGAIN.
      ALG = ALG - .1D0
      XN = N + IQ
      XTE = N + NMA
      E = E*XN/MIN(ALG*DBLE(N),XTE)
C     IF (MP.NE.0 .AND. LIN.GE.3) WRITE (NADV,99982) E, ALG
      GO TO 1660
C
C     CONSTRUCT  NEW  MESH  ***
C
 1700 J = 2
      I1 = 0
      DO 1740 I = 1, N
         DO 1720 JJ = 1, M
            I1 = I1 + 1
            UU(I1) = Y(JJ,I)
 1720    CONTINUE
 1740 CONTINUE
      SK(1) = A
      NPU = 2*NU + 3
      DO 1860 I = 1, N1
         KII1 = I*M
         IFI = J + IQJ(I)
         HI = HX(I)/DBLE(IQJ(I)+1)
         DO 1840 L = J, IFI
            SK(L) = X(I) + HI*DBLE(L-J+1)
            IF (IQJ(I).GT.0) GO TO 1780
            DO 1760 K = 1, M
               KII1PK = KII1 + K
               Y(K,L) = UU(KII1PK)
 1760       CONTINUE
            GO TO 1840
 1780       NP = NU + 2
            IF (I.GT.N-(NU+3)) NP = I - N + NPU
            IF (I.LE.NU+2) NP = I
            ITEMP = I
            CALL D02RAW(ITEMP,NPU,NP,AA(1,3),AA(1,2),X,NMAX,SK(L)
     *                  ,AA(1,5))
            DO 1820 K = 1, M
               YKI = ZERO
               DO 1800 JI = 1, NPU
                  IMNETC = M*(I-NP+JI-1) + K
                  YKI = YKI + AA(JI,3)*UU(IMNETC)
 1800          CONTINUE
               Y(K,L) = YKI
 1820       CONTINUE
 1840    CONTINUE
         J = IFI + 1
 1860 CONTINUE
C
      CASI = .FALSE.
      N = N + IQ
      N1 = N - 1
      DO 1880 I = 1, M
         MPNMPI = MPNM + I
         Y(I,N) = UU(MPNMPI)
 1880 CONTINUE
      H = ZERO
      DO 1900 I = 2, N1
         I1 = I - 1
         HX(I1) = SK(I) - SK(I1)
         IF (HX(I1).GT.H) H = HX(I1)
         X(I) = SK(I)
 1900 CONTINUE
C
      X(N) = B
      HX(N1) = X(N) - X(N1)
      IF (HX(N1).GT.H) H = HX(N1)
      IF (ABS(EPBAR-TEM*BMA).LT.1.D-5) EPBAR = 1.D10
      IF (ICON.GE.5) ICON = 15
      GO TO 260
C
C     *********** END OF MESH SELECTION *************
C     UPON ENTERING PASVA2 NOLD IS SET TO 1,
C     AND IT IS CHANGED TO N PRIOR TO ENTRANCE TO
C     THE SECOND MESH CHANGE. UPON ENTRY TO A
C     MESH CHANGE (AFTER THE VERY FIRST TIME)
C     NOLD = N . THUS, IF UPON EXIT OF THE MESH
C     SELECTION PROCEDURE STIL  NOLD=N THIS INDI-
C     CATES THAT NO CHANGE HAS TAKEN PLACE. THE
C     PARAMETER  ALG IS EQUAL TO 1.5 WHEN BEGI-
C     NING A MESH SELECTION AND WILL BE  =1.5  IF
C     AN ATTEMPT TO PUT MORE THAN  XNMA  POINTS
C     IN ONE STEP OF MESH REFINEMENT OCCURS.
C     NU  IS THE CORRECTION COUNTER, AND  NMAX IS
C     THE MAX. NUMBER OF MESH POINTS ALLOWED.
C     THE THREE NEXT LOGICAL INSTRUCTIONS DECIDE
C     WHICH PATH WILL BE TAKEN  ACCORDING TO THE
C     CONDITION OF THE VARIOUS PARAMETERS.
C     LABEL  1840  STARTS UNCONDITIONAL BISECTION
C     OF THE NET.  LABEL  1100  IS THE HEAD OF
C     THE CORRECTION PROCEDURE.
C     LABEL  1820    IS A CONDITIONAL BISECTION OF
C     THE NET (IF NU=0). SEQUENTIAL CONTINUATION
C     OCCURS WHEN  N=NOLD  .AND.  ALG=1.5
C     IT PROVIDES REENTRY TO THE MESH SELECTION
C     ALGORITHM WITH A LOWER LEVEL  E, IN AN
C     ATTEMPT TO FORCE POINTS INTO THE MESH.
C     IF THIS FAILS AND  NU .GT. 0, THEN WE RE-
C     ENTER TO THE MESH SELECTION ALGORITHM WITH
C     A METHOD OF ORDER  NU-1 (LABEL 1880).
C     ***********************************************
C
 1920 IF (NOLD.EQ.1 .AND. ALG.LT.1.5D0 .AND. 2*N-1.LE.NMAX) GO TO 1960
      IF (N.GT.NOLD) GO TO 1180
      IF (ALG.LT.1.45D0) GO TO 1940
      ALG = 1.4D0
      XN = N + IQ
      XXN = NMAX
      E = E*XN/MIN(ALG*DBLE(N),XXN)
C     IF (MP.NE.0 .AND. LIN.GE.3) WRITE (NADV,99982) E, ALG
      GO TO 1660
 1940 IF (NU.GT.0) GO TO 2000
C
C     BISECTION *****
C
      IF (MP.NE.0) THEN
         CALL X04BAF(NADV,' ')
         WRITE (REC,FMT=99986)
         CALL X04BAF(NADV,REC(1))
      END IF
      NTEM = 2*N - 1
      IF (NTEM.LE.NMAX) GO TO 1960
C
C     ... TOO MANY GRID POINTS  ...
C
      IFLAG = 2
      IF (LP.NE.0) THEN
         CALL X04BAF(NADV,' ')
         WRITE (REC,FMT=99988)
         CALL X04BAF(NADV,REC(1))
      END IF
      RETURN
 1960 DO 1980 I = 1, N1
         IQJ(I) = 1
 1980 CONTINUE
      IQ = N1
      ICON = 0
      GO TO 1700
 2000 NU = NU - 1
      ICON = 0
C
C     SINCE WE ARE GOING BACK WE MUST RECOMPUTE THE
C     ESTIMATE OF THE ERROR.
C
      MPNM = M*N1
      DO 2020 I = 1, MPNM
         SK(I) = RES(I) + SK(I)
 2020 CONTINUE
      IF (NU.GT.0) GO TO 2080
      DO 2060 I = 1, N1
         KI = (I-1)*M
         DO 2040 J = 1, M
            KIJ = KI + J
            RES(KIJ) = -SK(KIJ)
 2040    CONTINUE
 2060 CONTINUE
      GO TO 1600
 2080 CALL D02RAX(NU,2,2,N1,M,AA,X,NMAX,F,RES,MTNMAX,AA(1,4),AA(1,5)
     *            ,IERROR)
      DO 2120 I = 1, N1
         HI = HX(I)
         KI = (I-1)*M
         DO 2100 J = 1, M
            KIJ = KI + J
            RES(KIJ) = HI*RES(KIJ) - SK(KIJ)
 2100    CONTINUE
 2120 CONTINUE
      GO TO 1600
C     **** END OF MESH VARIATION ****
C
C     9982 FORMAT (10H   LEVEL =, 1P,E10.2, 8H   ALG =, 1P,E10.2)
C     9981 FORMAT (12H   N LEVEL =, 1P,E10.2, 16H   LOCAL ERROR =,
C     * 1PE10.2, 18H   MAX. SOLUTION =, 1PE10.2, 10H   EPBAR =,
C     * 1PE10.2)
99999 FORMAT (' NEWTON ITERATION HAS FAILED TO CONVERGE',/'  TRY A BET',
     *  'TER INITIAL APPROXIMATION AND/OR A FINER MESH')
99998 FORMAT ('  CONSIDER USING CONTINUATION')
99997 FORMAT ('  FOR MORE INFORMATION ABOUT FAILURE, SET IFAIL TO OBTA',
     *  'IN',/'  ADVISORY OUTPUT')
99996 FORMAT (' ZERO DERIVATIVES WITH RESPECT TO EPSILON',
     *  /'  EPSILON =',1P,D10.2,'   DELEPS =',1P,D10.2)
99995 FORMAT ('  CONTINUATION PARAMETER EPSILON =',1P,D10.2,'   DELEPS',
     *  ' =',1P,D10.2)
99994 FORMAT (' CONTINUATION ABANDONED --- DELEPS =',1P,D10.2)
99993 FORMAT ('  CORRECTION NUMBER ',I4,'   ESTIMATED MAXIMUM ERROR =',
     *  1P,D10.2)
99992 FORMAT ('  ESTIMATED ERROR BY COMPONENTS')
99991 FORMAT (2X,1P,7D10.2)
99990 FORMAT ('  NUMBER OF POINTS IN CURRENT MESH =',I5)
99989 FORMAT (' NEWTON ITERATION HAS REACHED ROUND-OFF LEVEL',/'  IF D',
     *  'ESIRED ACCURACY HAS NOT BEEN REACHED, THEN TOL IS TOO',/'  SM',
     *  'ALL FOR THIS PROBLEM AND THIS MACHINE PRECISION')
99988 FORMAT (' MORE THAN MNP MESH POINTS REQUIRED')
99987 FORMAT ('  MESH SELECTION')
99986 FORMAT ('  MESH BISECTION')
99985 FORMAT (' D02RAF MONITORING INFORMATION')
99984 FORMAT (' D02GBF MONITORING INFORMATION')
99983 FORMAT (' D02GAF MONITORING INFORMATION')
99982 FORMAT ('   NUMBER OF NEW POINTS ',I3)
99981 FORMAT ('   SINGULAR JACOBIAN')
      END
