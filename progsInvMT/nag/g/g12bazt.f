      SUBROUTINE G12BAZ(M,N,NS,Z,LDZ,AVG,T,ISZ,IP,ISI,IC,OFFL,OMEGA,DEV,
     *                  B,SE,SC,COV,SUR,NDMAX,IORD,NP,A,SI1,SS,SI2,WAK,
     *                  RES,ND,MAXIT,IPRINT,TOL,ITER,IERROR)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DEV, TOL
      INTEGER           IERROR, IP, IPRINT, ITER, LDZ, M, MAXIT, N, ND,
     *                  NDMAX, NS
      LOGICAL           OFFL
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IP), AVG(IP), B(IP), COV(IP*(IP+1)), OMEGA(*),
     *                  RES(N), SC(IP), SE(IP), SI1(IP), SI2(IP*(IP+1)),
     *                  SS(IP), SUR(NDMAX,*), T(N), WAK(N), Z(LDZ,M)
      INTEGER           IC(N), IORD(N), ISI(*), ISZ(M), NP(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  DEV0, DXI, S, SI, STEP, T1, T2, UFLOW, XMI
      INTEGER           I, IDI, II, INFO, IPOS, IPOS1, IPP, ISTEP,
     *                  ISTRAT, J, K, L, LL, MI, NLOOP, NOUT
      LOGICAL           CONVER
      CHARACTER*80      REC
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      EXTERNAL          X02AMF
C     .. External Subroutines ..
      EXTERNAL          DPPTRF, DPPTRI, DPPTRS, F06FBF, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, EXP, LOG, MOD, SQRT
C     .. Executable Statements ..
      IF (IPRINT.GT.0) CALL X04ABF(0,NOUT)
      IERROR = 0
      UFLOW = LOG(X02AMF())
      IPP = IP*(IP+1)/2
      IF (NS.EQ.0) THEN
         NLOOP = 1
         ISTRAT = 1
      ELSE
         NLOOP = NS
      END IF
C
C    Start Iterations
C
      CONVER = .FALSE.
      DEV0 = 0.0D0
      ITER = 0
   20 CONTINUE
      STEP = 1.0D0
      ISTEP = 0
   40 CONTINUE
      DEV = 0.0D0
      DO 80 I = 1, N
         IDI = IORD(I)
         IF (NS.GT.0) THEN
            ISTRAT = ISI(IDI)
         END IF
         IF (ISTRAT.GT.0) THEN
            IF (OFFL) THEN
               WAK(I) = OMEGA(IDI)
            ELSE
               WAK(I) = 0.0D0
            END IF
            K = 0
            DO 60 J = 1, M
               IF (ISZ(J).NE.0) THEN
                  K = K + 1
                  WAK(I) = WAK(I) + B(K)*(Z(IDI,J)-AVG(K))
               END IF
   60       CONTINUE
            IF (WAK(I).LT.UFLOW) THEN
               WAK(I) = 0.0D0
            ELSE IF (WAK(I).GT.-UFLOW) THEN
               IERROR = 4
               GO TO 580
            ELSE
               WAK(I) = EXP(WAK(I))
            END IF
         END IF
   80 CONTINUE
      CALL F06FBF(IPP,0.0D0,COV,1)
      CALL F06FBF(IP,0.0D0,SC,1)
      DO 280 LL = 1, NLOOP
         SI = 0.0D0
         CALL F06FBF(IPP,0.0D0,SI2,1)
         CALL F06FBF(IP,0.0D0,SI1,1)
         IPOS = N + 1
         DO 260 II = ND, 1, -1
            MI = 0
            CALL F06FBF(IP,0.0D0,SS,1)
            IPOS1 = IPOS - 1
            IPOS = NP(II)
            S = 0.0D0
            DO 180 I = IPOS, IPOS1
               IDI = IORD(I)
               IF (NS.GT.0) THEN
                  ISTRAT = ISI(IDI)
               END IF
               IF (ISTRAT.EQ.LL) THEN
                  K = 0
                  DO 100 J = 1, M
                     IF (ISZ(J).NE.0) THEN
                        K = K + 1
                        A(K) = Z(IDI,J) - AVG(K)
                     END IF
  100             CONTINUE
                  DXI = WAK(I)
                  SI = SI + DXI
                  L = 0
                  DO 140 K = 1, IP
                     SI1(K) = SI1(K) + DXI*A(K)
                     DO 120 J = 1, K
                        L = L + 1
                        SI2(L) = SI2(L) + DXI*A(K)*A(J)
  120                CONTINUE
  140             CONTINUE
                  IF (IC(IDI).EQ.0) THEN
                     IF (OFFL) S = S + OMEGA(IDI)
                     MI = MI + 1
                     DO 160 J = 1, IP
                        SS(J) = SS(J) + A(J)
  160                CONTINUE
                  END IF
               END IF
  180       CONTINUE
            IF (MI.NE.0 .AND. SI.GT.0.0D0) THEN
               DO 200 J = 1, IP
                  S = S + SS(J)*B(J)
  200          CONTINUE
               XMI = MI
               DEV = DEV + S - XMI*LOG(SI)
               L = 0
               DO 240 K = 1, IP
                  SC(K) = SC(K) - SS(K) + XMI*SI1(K)/SI
                  DO 220 J = 1, K
                     L = L + 1
                     T1 = SI2(L)/SI
                     T2 = (SI1(K)/SI)*(SI1(J)/SI)
                     COV(L) = COV(L) + XMI*(T1-T2)
  220             CONTINUE
  240          CONTINUE
            END IF
  260    CONTINUE
  280 CONTINUE
C
C         Check for convergance
C
      DEV = -2.0D0*DEV
      IF (ITER.GT.0) THEN
         IF (DEV.GT.DEV0 .AND. ISTEP.EQ.10) THEN
            IERROR = 6
            CONVER = .TRUE.
         ELSE IF (DEV.GT.DEV0) THEN
            ISTEP = ISTEP + 1
            STEP = STEP*0.5D0
            DO 300 I = 1, IP
               B(I) = B(I) - STEP*SE(I)
  300       CONTINUE
            GO TO 40
         ELSE IF ((DEV0-DEV).LT.(1.0D0+ABS(DEV0))*TOL) THEN
            CONVER = .TRUE.
         END IF
      END IF
      DEV0 = DEV
C
C       Iteration monitoring
C
      IF (IPRINT.GT.0) THEN
         IF (MOD(ITER,IPRINT).EQ.0) THEN
            WRITE (REC,FMT=99999) ITER, DEV
            CALL X04BAF(NOUT,REC)
            DO 320 I = 1, IP
               WRITE (REC,FMT=99998) I, B(I)
               CALL X04BAF(NOUT,REC)
  320       CONTINUE
         END IF
      END IF
C
C     Solve Newton-Raphson Equations
C
      CALL DPPTRF('U',IP,COV,INFO)
      IF (INFO.GT.0) THEN
         IERROR = 3
         GO TO 580
      END IF
      IF ( .NOT. CONVER) THEN
         ITER = ITER + 1
         IF (ITER.LT.MAXIT) THEN
            DO 340 J = 1, IP
               SE(J) = -SC(J)
  340       CONTINUE
            CALL DPPTRS('UPPER',IP,1,COV,SE,IP,INFO)
            DO 360 I = 1, IP
               B(I) = B(I) + SE(I)
  360       CONTINUE
            GO TO 20
         ELSE
            IERROR = 5
         END IF
      END IF
C
C          Compute variance-covariance matrix
C
      CALL DPPTRI('UPPER',IP,COV,INFO)
      IF (INFO.NE.0) THEN
         IERROR = 3
         GO TO 580
      ELSE
C
C           Compute S.E.
C
         DO 380 K = 1, IP
            SE(K) = SQRT(COV(K*(K+1)/2))
  380    CONTINUE
      END IF
C
C          Compute hazard functions
C
      DO 420 J = 1, NLOOP
         DO 400 I = 1, ND
            SUR(I,J) = 0.0D0
  400    CONTINUE
  420 CONTINUE
      DO 440 I = 1, NP(1) - 1
         IDI = IORD(I)
         RES(IDI) = 0.0D0
  440 CONTINUE
      DO 560 LL = 1, NLOOP
         IPOS = N + 1
         SI = 0.0D0
         DO 480 II = ND, 1, -1
            MI = 0
            IPOS1 = IPOS - 1
            IPOS = NP(II)
            DO 460 I = IPOS, IPOS1
               IDI = IORD(I)
               IF (NS.GT.0) THEN
                  ISTRAT = ISI(IDI)
               END IF
               IF (ISTRAT.EQ.LL) THEN
                  SI = SI + WAK(I)
                  IF (IC(IDI).EQ.0) THEN
                     MI = MI + 1
                  END IF
               END IF
  460       CONTINUE
            IF (MI.NE.0 .AND. SI.GT.0.0D0) THEN
               SUR(II,LL) = DBLE(MI)/SI
            END IF
  480    CONTINUE
         SI = 0.0D0
         DO 500 II = 1, ND
            SI = SI + SUR(II,LL)
            SUR(II,LL) = SI
  500    CONTINUE
         IPOS = N + 1
         DO 540 II = ND, 1, -1
            IPOS1 = IPOS - 1
            IPOS = NP(II)
            DO 520 I = IPOS, IPOS1
               IDI = IORD(I)
               IF (NS.GT.0) THEN
                  ISTRAT = ISI(IDI)
               END IF
               IF (ISTRAT.EQ.LL) THEN
                  RES(IDI) = WAK(I)*SUR(II,LL)
               END IF
  520       CONTINUE
            SUR(II,LL) = EXP(-SUR(II,LL))
  540    CONTINUE
  560 CONTINUE
  580 RETURN
C
99999 FORMAT (' Iteration ',I16,' Deviance = ',D13.5)
99998 FORMAT (' Beta(',I6,') = ',D13.5)
      END
