      SUBROUTINE D02SAX(P,N,A,B,U,V,D,IV,ERRP,PF,F,W,MONIT)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 9 REVISED. IER-311 (SEP 1981).
C     MARK 11 REVISED. IER-419 (FEB 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 15 REVISED. IER-903 (APR 1991).
C     NEWTON WITH SVD.
C     MONIT
C     .. Scalar Arguments ..
      INTEGER           IV, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IV,2), B(IV,2), D(N), ERRP(N), F(N), P(N),
     *                  PF(N), U(IV,N), V(IV,N), W(IV,6)
C     .. Subroutine Arguments ..
      EXTERNAL          MONIT
C     .. Scalars in Common ..
      DOUBLE PRECISION  COUT12, COUT13, DM, EPS, EPSFAC, MACHEP, PNORM,
     *                  PNORM1, SQEPS
      INTEGER           COUNT, IC, ICASE, IEPS, IFAIL1, IFLAG, ISTATE
C     .. Local Scalars ..
      DOUBLE PRECISION  DMIN, S
      INTEGER           I, J
C     .. External Subroutines ..
      EXTERNAL          D02SAT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, DBLE
C     .. Common blocks ..
      COMMON            /AD02SA/PNORM, PNORM1, EPS, MACHEP, SQEPS, DM,
     *                  EPSFAC, ISTATE, IFLAG, COUNT, IFAIL1, IEPS,
     *                  ICASE
      COMMON            /BD02SA/COUT12, COUT13, IC
C     .. Executable Statements ..
      IF (N.GT.0 .AND. IV.GE.N) GO TO 20
C     INPUT ERROR
      IFLAG = -1
      RETURN
   20 IF (ISTATE.NE.0) GO TO 40
      IFLAG = 1
      EPS = 0.D0
      EPSFAC = SQEPS**(1.0D0/4.0D0)
      IEPS = 0
      GO TO 80
   40 IF (ISTATE.LT.0) GO TO 640
      IF (ISTATE-6) 300, 80, 60
C     INPUT ERROR IN COMMON
   60 IFLAG = -1
      RETURN
C     ENTRY AFTER CALCULATING SVD OF JACOBIAN
   80 DM = 0.D0
      DMIN = D(1)
      DO 100 I = 1, N
         DM = MAX(DM,D(I))
         DMIN = MIN(DMIN,D(I))
         W(I,6) = D(I)
  100 CONTINUE
      IF (DMIN.LT.SQEPS*DM .AND. EPS.EQ.0.D0) EPS = SQEPS/EPSFAC
C     PERTURB THE SVD
  120 DO 140 I = 1, N
         D(I) = D(I)/(D(I)*D(I)+0.25D0*EPS*EPS*DM*DM)
  140 CONTINUE
      DO 180 I = 1, N
         S = 0.D0
         DO 160 J = 1, N
            S = S + U(J,I)*F(J)
  160    CONTINUE
         W(I,4) = S*D(I)
  180 CONTINUE
      IF (ICASE.EQ.3) CALL D02SAT(P,A,B,IV,0)
      PNORM = 0.D0
      DO 220 I = 1, N
         S = 0.D0
         DO 200 J = 1, N
C Changed to reflect mod in D02SAW - R.W.Brankin, NAG, Feb 20 1991
C           S = S + V(I,J)*W(J,4)
            S = S + V(J,I)*W(J,4)
  200    CONTINUE
         W(I,1) = S
         W(I,5) = F(I)
         S = S/MAX(ABS(P(I)),PF(I))
         PNORM = PNORM + S*S
  220 CONTINUE
      CALL MONIT(ISTATE,IFLAG,IFAIL1,P,N,F,PNORM,PNORM,EPS,W(1,6))
      IF (ISTATE.NE.0) GO TO 260
      DO 240 I = 1, N
         IF (ABS(W(I,1)).GT.ERRP(I)*MAX(ABS(P(I)),PF(I))) GO TO 260
  240 CONTINUE
      IFLAG = 0
      RETURN
  260 CONTINUE
      ISTATE = 1
      DO 280 I = 1, N
         W(I,2) = P(I)
         P(I) = P(I) - W(I,1)
  280 CONTINUE
      IF (ICASE.EQ.3) CALL D02SAT(P,A,B,IV,1)
      RETURN
C     ENTRY AFTER RECALCULATING RESIDUAL
  300 DO 340 I = 1, N
         S = 0.D0
         DO 320 J = 1, N
            S = S + U(J,I)*F(J)
  320    CONTINUE
         W(I,4) = S*D(I)
  340 CONTINUE
      IF (ICASE.EQ.3) CALL D02SAT(P,A,B,IV,0)
      PNORM1 = 0.D0
      DO 380 I = 1, N
         S = 0.D0
         DO 360 J = 1, N
C Changed to reflect mod in D02SAW - R.W.Brankin, NAG, Feb 20 1991
C           S = S + V(I,J)*W(J,4)
            S = S + V(J,I)*W(J,4)
  360    CONTINUE
         W(I,3) = S
         S = S/MAX(ABS(P(I)),PF(I))
         PNORM1 = PNORM1 + S*S
  380 CONTINUE
      IF (ISTATE.EQ.1) COUT13 = COUT12
      CALL MONIT(ISTATE,IFLAG,IFAIL1,P,N,F,PNORM,PNORM1,EPS,W(1,6))
      IF (PNORM1.GE.0.9D0*PNORM) GO TO 500
      IF (ISTATE.GT.1) GO TO 480
      IF (EPS.GE.SQEPS) GO TO 460
C     TAKES FULL NEWTON STEP
      DO 400 I = 1, N
         IF (ABS(W(I,1)).GT.ERRP(I)*MAX(ABS(P(I)),PF(I))) GO TO 420
  400 CONTINUE
C     CONVERGENCE CRITERION SATISFIED
      IFLAG = 0
      RETURN
  420 IF (PNORM1.GE.0.0625D0*(1.D0-1.D0/DBLE(N+1))*PNORM) GO TO 460
C     NO RECALCULATION OF SVD
      PNORM = PNORM1
      IFLAG = 2
      DO 440 I = 1, N
         W(I,1) = W(I,3)
         W(I,2) = P(I)
         W(I,5) = F(I)
         P(I) = P(I) - W(I,1)
  440 CONTINUE
      IF (ICASE.EQ.3) CALL D02SAT(P,A,B,IV,1)
      RETURN
C     OUT OF TROUBLE
  460 EPS = 0.D0
      IEPS = 0
C     RE-EVALUATE SVD
  480 IFLAG = 1
      ISTATE = 6
      RETURN
C     HALVE THE NEWTON STEP
  500 IF (ISTATE.EQ.5) GO TO 560
  520 ISTATE = ISTATE + 1
      IF (IFLAG.EQ.2 .AND. ISTATE.EQ.3) GO TO 560
      DO 540 I = 1, N
         W(I,1) = 0.5D0*W(I,1)
         P(I) = W(I,2) - W(I,1)
  540 CONTINUE
      IF (ICASE.EQ.3) CALL D02SAT(P,A,B,IV,1)
      RETURN
C     REJECT NEWTON STEP
  560 DO 580 I = 1, N
         F(I) = W(I,5)
         P(I) = W(I,2)
  580 CONTINUE
      IF (ICASE.EQ.3) CALL D02SAT(P,A,B,IV,1)
      ISTATE = 6
      IF (IFLAG.EQ.1) GO TO 600
C     RE-EVALUATE SVD
      IFLAG = 1
      COUT12 = COUT13
      RETURN
C     CHANGE SINGULAR VALUES
  600 EPS = EPS/EPSFAC
      IF (EPS.EQ.0.D0) EPS = SQEPS/EPSFAC
      IEPS = IEPS + 1
      IF (IEPS.EQ.4) GO TO 660
      DO 620 I = 1, N
         D(I) = W(I,6)
  620 CONTINUE
      GO TO 120
C     ERROR ENTRY
  640 IF (ICASE.EQ.3) CALL D02SAT(P,A,B,IV,0)
      CALL MONIT(ISTATE,IFLAG,IFAIL1,P,N,F,PNORM,0.D0,EPS,W(1,6))
      ISTATE = -ISTATE
      IF (ISTATE.EQ.7) RETURN
      IF (ISTATE-5) 520, 560, 600
C     NO CONVERGENCE
  660 IFLAG = -2
      RETURN
      END
