      SUBROUTINE G13DCF(K,N,P,Q,MEAN,X,N4,QQ,IK,W,PARHLD,CONDS,IPRINT,
     *                  CGETOL,MAXCAL,ISHOW,NITER,LOGL,V,G,DISP,IDISP,
     *                  W2,LW,IW,LIW,IFAIL)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1111 (JUL 1993).
C     MARK 17 REVISED. IER-1679 (JUL 1995).
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13DCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CGETOL, LOGL
      INTEGER           IDISP, IFAIL, IK, IPRINT, ISHOW, K, LIW, LW,
     *                  MAXCAL, N, N4, NITER, P, Q
      LOGICAL           CONDS, MEAN
C     .. Array Arguments ..
      DOUBLE PRECISION  DISP(IDISP,N4), G(N4), QQ(IK,K), V(IK,N),
     *                  W(IK,N), W2(LW), X(N4)
      INTEGER           IW(LIW)
      LOGICAL           PARHLD(N4)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ADDLOG, COND, DETQQ, EXPP, FBIG, LMAX, NORM,
     *                  SMALL, SSQ, SUMSSQ, XTOL
      INTEGER           IFAILX, INTEG, IPP, IQQ, ITN, K3, K4, K5, K6,
     *                  KK, KMAT, KR, KZ, LEW6, LEW7, LIW1, LIWW, LP,
     *                  LW1, LW10, LW11, LW12, LW13, LW14, LW15, LW16,
     *                  LW17, LW18, LW19, LW2, LW20, LW21, LW3, LW4,
     *                  LW5, LW6, LW7, LW8, LW9, LWW, NN, NNN, R
      LOGICAL           ACONDS, AMEAN, ANOTT, FULLP, FULLQ, NOPRIN
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, EPSA, ETA, F, FEST, STEPMX, SUM, XXTOL
      INTEGER           I, I2, I3, I6, IBOUND, IDIAG, IFAILY, IFLAG,
     *                  IFP, INTYPE, IPR, ITT, IWARN, J, K7, L, LAGNUM,
     *                  LEW1, LEW10, LEW11, LEW2, LEW3, LEW4, LEW5,
     *                  LEW8, LEW9, LH, LIW2, MODE, MSGLVL, N2
      LOGICAL           DELTA, LOCSCH, SETQQ
C     .. Local Arrays ..
      CHARACTER*80      P01REC(5), REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF, X02AJF, X02ALF, X02AMF
      INTEGER           P01ABF
      EXTERNAL          X01AAF, X02AJF, X02ALF, X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          E04HBZ, E04JBL, E04JBQ, E04XAF, F01AAZ, G13DCP,
     *                  G13DCQ, G13DCR, G13DCS, G13DCT, G13DCU, G13DCW,
     *                  G13DCX, G13DCZ, X04ABF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, LOG, MAX, MIN, DBLE, SQRT
C     .. Common blocks ..
      COMMON            /AG13DC/LW6, LW7, LW9, LW10, LW15, LW19, LW20,
     *                  LW21, LIW1
      COMMON            /BG13DC/ADDLOG, LMAX, COND, NORM, KK, IPP, IQQ,
     *                  K4, K5, LW14, LW18, LP, K6, ITN, LEW6, LEW7,
     *                  LW2, LW1, LW3, LW4, LW5, LW8, LW11, LW12, LW13,
     *                  LW16, LW17, KR, AMEAN, NOPRIN, FULLP, FULLQ
      COMMON            /CG13DC/XTOL, SSQ, SUMSSQ, EXPP, DETQQ, FBIG,
     *                  SMALL, NN, R, INTEG, K3, KZ, KMAT, IFAILX,
     *                  ANOTT, ACONDS
      COMMON            /DG13DC/NNN, LIWW, LWW
C     .. Executable Statements ..
C
C     INITIALISE SCALARS IN COMMON BLOCKS
C
      CALL X04ABF(0,LP)
      IF (CGETOL.LT.X02AJF()) THEN
         XXTOL = 10.0D0*SQRT(X02AJF())
         WRITE (REC,FMT=99997) CGETOL, XXTOL
         CALL X04BAY(LP,2,REC)
      ELSE
         XXTOL = CGETOL
      END IF
      CONDS = .NOT. CONDS
      ACONDS = CONDS
      AMEAN = MEAN
      IPP = P
      IQQ = Q
      NN = N
      KK = K
      LIWW = LIW
      LWW = LW
      ADDLOG = (DBLE(N)*DBLE(K)/2.0D0)*LOG(2.0D0*X01AAF(ADDLOG))
      SMALL = 10.0D0*X02AMF()
C
      IF (CONDS) THEN
         XTOL = X02ALF()/100.0D0
      ELSE
         XTOL = 0.001D0
      END IF
C
      R = MAX(P,Q)
      KR = K*R
C
      K3 = K*K
      K4 = P*K3
      K5 = K4 + Q*K3
      K6 = K5
      IF (MEAN) K6 = K5 + K
C
C     TEST FOR ERRORS IN THE INPUT DATA
C
      IFAILX = 0
      IF (K.LT.1) IFAILX = 1
      IF (N4.LT.1) IFAILX = 1
      IF (N*K.LE.N4+K*(K+1)/2) IFAILX = 1
      IF (P.LT.0) IFAILX = 1
      IF (Q.LT.0) IFAILX = 1
      IF ((P.EQ.0) .AND. (Q.EQ.0)) IFAILX = 1
      IF (N4.NE.K6) IFAILX = 1
      IF (IK.LT.K) IFAILX = 1
      IF (MAXCAL.LT.1) IFAILX = 1
      IF ((ISHOW.LT.0) .OR. (ISHOW.GT.2)) IFAILX = 1
      IF (IDISP.LT.N4) IFAILX = 1
C
C     SET UP PARTITIONING OF WORKSPACE ARRAYS IW AND W2
C
      LIW1 = 3
      KMAT = K*K*(R+1)
      KZ = K*K*(R+1)
C
      N2 = N4 + K*(K+1)/2
      LW1 = 9*N2 + 1
      LW2 = LW1 + K*(P*K+1)
      LW3 = LW2 + K*(Q*K+1)
      LW4 = LW3 + K*KR
      LW5 = LW4 + K*KR
      LW6 = LW5 + K*K
      LW7 = LW6 + K*N
      LW8 = LW7 + K*N
      LW9 = LW8 + K*K
      LW10 = LW9 + KZ
      LW11 = LW10 + KMAT*KMAT
      LW12 = LW11 + K*(K+1)
      LW13 = LW12 + K*K
      LW14 = LW13 + K*KR
      LW15 = LW14 + K*K
      LW16 = LW15 + K*(K+1)
      LW17 = LW16 + K*K*(P+1)
      LW18 = LW17 + K*K
      LW19 = LW18 + K
      LW20 = LW19 + KR
      LW21 = LW20 + KR
C
C     SET UP PARTITIONING OF THE E04 ARRAY ARGUMENTS
C
      LEW1 = LW21 + K*N
      LEW2 = LEW1 + N2
      LEW3 = LEW2 + N2
      LEW4 = LEW3 + N2
      LH = N2*(N2-1)/2
      LEW5 = LEW4 + LH
      LEW6 = LEW5 + N2
      LEW7 = LEW6 + KR*KR
      LEW8 = LEW7 + KR*KR
      LEW9 = LEW8 + N2*N2
      LEW10 = LEW9 + N2*N2
      LEW11 = LEW10 + N2
      LIW2 = LIW1 + KR
      IF (LW.LT.LEW11+N2-1) THEN
         IFAILX = 1
         WRITE (P01REC,FMT=99996) LW, LEW11 + N2 - 1
         IFAIL = P01ABF(IFAIL,IFAILX,SRNAME,1,P01REC)
         RETURN
      END IF
C
      IF (LIW.LT.LIW2+N2-1) THEN
         IFAILX = 1
         WRITE (P01REC,FMT=99995) LIW, LIW2 + N2 - 1
         IFAIL = P01ABF(IFAIL,IFAILX,SRNAME,1,P01REC)
         RETURN
      END IF
C
      IF (IFAILX.NE.0) THEN
         WRITE (P01REC,FMT=99991) K, N, P, Q, N4, IK, MAXCAL, ISHOW,
     *     IDISP, MEAN
         IFAIL = P01ABF(IFAIL,IFAILX,SRNAME,5,P01REC)
         RETURN
      END IF
C
C     SET IBOUND
C
      IBOUND = 1
      LAGNUM = 0
      DO 20 I = 1, N4
         IF (PARHLD(I)) THEN
            IBOUND = 0
            LAGNUM = LAGNUM + 1
         END IF
   20 CONTINUE
C
      IF (IBOUND.EQ.1) THEN
         FULLP = .TRUE.
      ELSE
         FULLP = .TRUE.
         IF (P.GT.0) THEN
            DO 40 I = 1, P*K*K
               IF (PARHLD(I)) FULLP = .FALSE.
   40       CONTINUE
         END IF
      END IF
C
      IF (IBOUND.EQ.1) THEN
         FULLQ = .TRUE.
      ELSE
         FULLQ = .TRUE.
         IF (Q.GT.0) THEN
            DO 60 I = K4 + 1, K5
               IF (PARHLD(I)) FULLQ = .FALSE.
   60       CONTINUE
         END IF
      END IF
C
C     TEST WHETHER SIGMA IS POSITIVE DEFINITE AND EXIT IF IT ISN'T
C
C     FIRST TEST WHETHER ALL ELEMENTS OF SIGMA ARE ZERO, IN WHICH CASE
C     WE NEED TO INITIALISE QQ TO THE SAMPLE COVARIANCE MATRIX
C
      SETQQ = .TRUE.
      DO 100 J = 1, K
         DO 80 I = J, K
            IF (QQ(I,J).NE.0.0D0) SETQQ = .FALSE.
   80    CONTINUE
  100 CONTINUE
C
      IF (SETQQ) CALL G13DCP(QQ,IK,K,W,N,W2)
C
      CALL G13DCW(QQ,IK,K,DELTA,V)
      IF ( .NOT. DELTA) THEN
         IFAILX = 2
         WRITE (P01REC,FMT=99994)
     *     '** INITIAL ESTIMATE OF SIGMA IS NOT POSITIVE-DEFINITE'
         IFAIL = P01ABF(IFAIL,IFAILX,SRNAME,1,P01REC)
         RETURN
      END IF
C
C     SCALE THE CHOLESKI DECOMPOSITION OF SIGMA
C
      LMAX = ABS(V(1,1))
C
      IF (K.EQ.1) GO TO 160
      DO 140 I = 2, K
         DO 120 J = 1, I
            IF (LMAX.LT.ABS(V(I,J))) LMAX = ABS(V(I,J))
  120    CONTINUE
  140 CONTINUE
C
  160 DO 200 I = 1, K
         DO 180 J = 1, I
            W2(LEW10-1+K6+(I-1)*I/2+J) = V(I,J)/LMAX
  180    CONTINUE
  200 CONTINUE
C
C     SET V TO THE IDENTITY MATRIX
C
      DO 240 I = 1, K
         DO 220 J = 1, K
            V(I,J) = 0.0D0
            IF (I.EQ.J) V(I,J) = 1.0D0
  220    CONTINUE
  240 CONTINUE
C
C     TEST WHETHER THE STARTING POINT IS INSIDE THE STATIONARITY
C     REGION
C
      IF (P.GT.0) THEN
         DO 300 L = 1, P
            DO 280 I = 1, K
               DO 260 J = 1, K
                  W2(LW1-1+(L-1)*K*K+(J-1)*K+I) = X((L-1)*K*K+(I-1)*K+J)
  260          CONTINUE
  280       CONTINUE
  300    CONTINUE
         CALL G13DCX(P,K,W2(LW1),W2(LW10),KMAT,R,W2(LW19),W2(LW20),
     *               IW(LIW1),KR,DELTA)
         IF ( .NOT. DELTA) THEN
            IFAILX = 2
            WRITE (P01REC,FMT=99994)
     *   '** INITIAL AR PARAMETER ESTIMATES OUTSIDE STATIONARITY REGION'
            IFAIL = P01ABF(IFAIL,IFAILX,SRNAME,1,P01REC)
            RETURN
         END IF
      END IF
C
C      TEST WHETHER THE STARTING POINT IS INSIDE THE INVERTIBILITY
C      REGION
C
      IF (Q.GT.0) THEN
         DO 360 L = 1, Q
            DO 340 I = 1, K
               DO 320 J = 1, K
                  W2(LW2-1+(L-1)*K*K+(J-1)*K+I) = X(P*K*K+(L-1)
     *              *K*K+(I-1)*K+J)
  320          CONTINUE
  340       CONTINUE
  360    CONTINUE
         CALL G13DCX(Q,K,W2(LW2),W2(LW10),KMAT,R,W2(LW19),W2(LW20),
     *               IW(LIW1),KR,DELTA)
         IF ( .NOT. DELTA) THEN
            IFAILX = 2
            WRITE (P01REC,FMT=99994)
     *  '** INITIAL MA PARAMETER ESTIMATES OUTSIDE INVERTIBILITY REGION'
            IFAIL = P01ABF(IFAIL,IFAILX,SRNAME,1,P01REC)
            RETURN
         END IF
      END IF
C
      IF (P.GT.0) THEN
C
C        IF WE ARE REPARAMETERISING ON THE AR SIDE OF THE ARMA
C        EQUATION FIND THE A(J)'S
C
         IF (FULLP) THEN
            I6 = 0
            CALL G13DCU(I6,P,X,N4,K,KR,KMAT,W2(LW1),V,IK,W2(LW3),
     *                  W2(LW10),W2(LW9),KZ,W2(LW16),W2(LEW6),W2(LEW7),
     *                  W2(LW5),W2(LW8),W2(LW11),W2(LW12),W2(LW14),
     *                  W2(LW15),W2(LEW10),IFAILX)
            IF (IFAILX.EQ.2) THEN
               WRITE (P01REC,FMT=99994)
     *   '** INITIAL AR PARAMETER ESTIMATES OUTSIDE STATIONARITY REGION'
               IFAIL = P01ABF(IFAIL,IFAILX,SRNAME,1,P01REC)
               RETURN
            END IF
         END IF
      END IF
C
      IF (Q.GT.0) THEN
C
C        IF WE ARE REPARAMETERISING ON THE MA SIDE OF THE ARMA
C        EQUATION FIND THE A(J)'S
C
         IF (FULLQ) THEN
            CALL G13DCU(K4,Q,X,N4,K,KR,KMAT,W2(LW2),V,IK,W2(LW3),
     *                  W2(LW10),W2(LW9),KZ,W2(LW13),W2(LEW6),W2(LEW7),
     *                  W2(LW17),W2(LW8),W2(LW11),W2(LW12),W2(LW5),
     *                  W2(LW15),W2(LEW10),IFAILX)
            IF (IFAILX.EQ.2) THEN
               WRITE (P01REC,FMT=99994)
     *  '** INITIAL MA PARAMETER ESTIMATES OUTSIDE INVERTIBILITY REGION'
               IFAIL = P01ABF(IFAIL,IFAILX,SRNAME,1,P01REC)
               RETURN
            END IF
         END IF
      END IF
C
C     SET UP THE PARAMETERS FOR E04JBL
C
      IFLAG = 0
      LOCSCH = .TRUE.
      INTYPE = 1
      FEST = 0.5D0
      ETA = 0.5D0
      IF (N2-LAGNUM.EQ.1) ETA = 0.0D0
      STEPMX = 100000.0D0
      ANOTT = (P.GT.Q)
C
C     PUT TIME SERIES COLUMN-WISE INTO W2
C
      DO 400 J = 1, N
         DO 380 I = 1, K
            W2(LW6-1+(J-1)*K+I) = W(I,J)
  380    CONTINUE
  400 CONTINUE
C
      IF ( .NOT. MEAN) GO TO 480
C
C     CALCULATE MEAN OF EACH SERIES AND STORE IN ARRAY MEANS
C
      DO 460 I = 1, K
         IF (PARHLD(K5+I)) THEN
            W2(LW18-1+I) = 0.0D0
            W2(LEW10-1+K5+I) = X(K5+I)
         ELSE
            IF (X(K5+I).EQ.0.0D0) THEN
               SUM = 0.0D0
               DO 420 J = 1, N
                  SUM = SUM + W(I,J)
  420          CONTINUE
               W2(LW18-1+I) = SUM/DBLE(N)
               W2(LEW10-1+K5+I) = 0.0D0
            ELSE
               SUM = 0.0D0
               DO 440 J = 1, N
                  SUM = SUM + W(I,J)
  440          CONTINUE
               W2(LW18-1+I) = SUM/DBLE(N)
               W2(LEW10-1+K5+I) = X(K5+I) - W2(LW18-1+I)
            END IF
         END IF
  460 CONTINUE
C
C     IF ANY PARAMETERS ARE TO BE HELD EQUAL TO ZERO
C     WE MUST NOTE WHICH ONES THEY ARE
C
  480 IF (IBOUND.EQ.1) GO TO 540
C
      DO 500 I = 1, N2
         W2(LEW1-1+I) = 1000000.0D0
         W2(LEW2-1+I) = -1000000.0D0
  500 CONTINUE
      DO 520 I = 1, N4
         IF (PARHLD(I)) W2(LEW1-1+I) = X(I)
         IF (PARHLD(I)) W2(LEW2-1+I) = X(I)
  520 CONTINUE
C
C     CALL G13DCR SO THAT WE MAY CALCULATE THE NEGATIVE OF THE
C     LOG LIKELIHOOD FUNCTION AT THE USER SUPPLIED STARTING POINT
C
  540 INTEG = 0
C
C     IF FULLP = FALSE (AND FULLQ =FALSE) PUT PHI (AND THETA)
C     PART OF X ARRAY INTO W2(LEW10)
C
      IF ((P.GT.0) .AND. ( .NOT. FULLP)) THEN
         DO 560 I = 1, K4
            W2(LEW10-1+I) = X(I)
  560    CONTINUE
      END IF
C
      IF ((Q.GT.0) .AND. ( .NOT. FULLQ)) THEN
         DO 580 I = K4 + 1, K5
            W2(LEW10-1+I) = X(I)
  580    CONTINUE
      END IF
C
      IFLAG = 0
C
      CALL G13DCR(IFLAG,N2,W2(LEW10),F,W2(LEW11),IW,LIW,W2,LW)
C
      IF ((IFAILX.EQ.2) .OR. (IFLAG.EQ.-1)) THEN
         IFAILX = 2
         WRITE (P01REC,FMT=99993)
     *     '** STARTING POINT IS TOO CLOSE TO BOUNDARY OF ',
     *     'ADMISSABILITY REGION'
         IFAIL = P01ABF(IFAIL,IFAILX,SRNAME,1,P01REC)
         RETURN
      END IF
C
      INTEG = 1
      IFAILY = 1
C
C     SET UP DIFFERENCE INTERVALS FOR CALCULATING DERIVATIVES
C
      CALL E04HBZ(N2,G13DCR,W2(LEW10),J,W2(LEW3),W2(LEW4),LH,W2(LEW5),F,
     *            W2(LEW11),IW,LIW,W2,LW,IFAILY)
C
      IF (IFAILY.NE.0) THEN
         IF (IFAILY.EQ.-1) IFAILY = 2
         IFAILY = IFAILY + 1
         WRITE (P01REC,FMT=99993)
     *     '** MINIMIZATION ROUTINE CANNOT COMPUTE AN ACCURATE ESTIMATE'
     *     , ' OF THE GRADIENT VECTOR'
         IFAIL = P01ABF(IFAIL,IFAILY,SRNAME,1,P01REC)
         RETURN
      END IF
C
      IFAILY = 1
C
      IPR = IPRINT
      NOPRIN = .FALSE.
      IF (IPRINT.LT.0) THEN
         NOPRIN = .TRUE.
         IPRINT = 0
      END IF
C
C     CALL E04JBL TO MAXIMISE THE LOG LIKELIHOOD FUNCTION
C
      CALL E04JBL(N2,G13DCR,G13DCQ,IPRINT,LOCSCH,INTYPE,E04JBQ,MAXCAL,
     *            ETA,XXTOL,STEPMX,FEST,W2(LEW3),IBOUND,W2(LEW2),
     *            W2(LEW1),W2(LEW10),W2(LEW4),LH,W2(LEW5),IW(LIW2),F,
     *            W2(LEW11),IW,LIW,W2,LW,IFAILY)
C
      IF (IFAILY.EQ.1) THEN
         IFAILY = IFAILY + 2
         IPRINT = IPR
         WRITE (P01REC,FMT=99993)
     *     '** MINIMIZATION ROUTINE CANNOT COMPUTE AN ACCURATE ',
     *     'ESTIMATE OF GRADIENT VECTOR'
         IFAIL = P01ABF(IFAIL,IFAILY,SRNAME,1,P01REC)
         RETURN
      END IF
C
C     TEST WHETHER E04JBL HAS FAILED TO RECOGNISE THE FINAL POINT
C     AS THE SOLUTION
C
      IF ((IFAILY.EQ.3) .OR. (IFAILY.EQ.5)) THEN
         EPS = 10.0D0*X02AJF()
         IF (NORM*NORM.LT.EPS) IFAILY = 0
      END IF
C
      SETQQ = (IFAILY.EQ.0)
      IPRINT = IPR
      NITER = ITN
C
C     COMPUTE THE 'FULL' LOG LIKELIHOOD FUNCTION AFTER EXIT FROM
C     E04JBL
C
      INTEG = 2
      EPSA = XTOL
      CALL G13DCR(IFLAG,N2,W2(LEW10),F,W2(LEW11),IW,LIW,W2,LW)
C
C     EXTRACT V ARRAY FROM W2
C
      DO 620 J = 1, N
         DO 600 I = 1, K
            V(I,J) = W2(LW7-1+(J-1)*K+I)
  600    CONTINUE
  620 CONTINUE
C
      LOGL = -F - ADDLOG
C
C     RECONSTRUCT SIGMA MATRIX
C
      DO 680 I = 1, K
         DO 660 J = 1, I
            QQ(I,J) = 0.0D0
            DO 640 I3 = 1, K
               IF (MIN(I,J).GE.I3) QQ(I,J) = QQ(I,J) +
     *             W2(LEW10-1+K6+(I-1)*I/2+I3)*W2(LEW10-1+K6+(J-1)
     *             *J/2+I3)
  640       CONTINUE
            QQ(I,J) = LMAX*LMAX*QQ(I,J)
            QQ(J,I) = QQ(I,J)
C
  660    CONTINUE
  680 CONTINUE
C
C     SET ELEMENTS OF X TO THEIR ML ESTIMATES
C
      DO 700 I = 1, K6
         X(I) = W2(LEW10-1+I)
  700 CONTINUE
C
      IF (FULLP .AND. (P.GT.0)) THEN
         DO 760 L = 1, P
            DO 740 J = 1, K
               DO 720 I = 1, K
                  X((L-1)*K*K+(I-1)*K+J) = W2(LW1-1+(L-1)*K*K+(J-1)*K+I)
  720          CONTINUE
  740       CONTINUE
  760    CONTINUE
         DO 780 I = 1, P*K*K
            W2(LEW10-1+I) = X(I)
  780    CONTINUE
      END IF
C
      IF (FULLQ .AND. (Q.GT.0)) THEN
         DO 840 L = 1, Q
            DO 820 J = 1, K
               DO 800 I = 1, K
                  X(P*K*K+(L-1)*K*K+(I-1)*K+J) = W2(LW2-1+(L-1)
     *              *K*K+(J-1)*K+I)
  800          CONTINUE
  820       CONTINUE
  840    CONTINUE
         DO 860 I = 1, Q*K*K
            W2(LEW10-1+P*K*K+I) = X(P*K*K+I)
  860    CONTINUE
      END IF
C
      IF ((IFAILY.NE.0) .AND. (IFAIL.LE.1)) THEN
         IFAILY = IFAILY + 2
         IF ((IFAILY.EQ.6) .OR. (IFAILY.EQ.7)) IFAILY = 5
         IF (IFAILY.EQ.4) WRITE (P01REC,FMT=99992) ' ** MAXCAL (',
     *       MAXCAL, ' ) LIKELIHOOD EVALUATIONS HAVE BEEN MADE'
         IF (IFAILY.EQ.5) WRITE (P01REC,FMT=99994)
     *       '** CONDITIONS FOR A SOLUTION HAVE NOT ALL BEEN MET'
         IFAIL = P01ABF(IFAIL,IFAILY,SRNAME,1,P01REC)
      END IF
C
C     COMPUTE HESSIAN MATRIX AT SOLUTION POINT
C
      K7 = N4
      IF (IBOUND.EQ.0) K7 = N4 - LAGNUM
C
      INTEG = 3
      FULLP = .FALSE.
      FULLQ = .FALSE.
      MODE = 2
      IDIAG = 1
      MSGLVL = 0
      EPSA = X02AJF()*(1.0D0+ABS(F))
C
      NNN = N4
C
C     COPY W2(LEW10) (REDUCED) ONTO W2(1) AND USE W2(K7+1) FOR
C     THE GRADIENT VECTOR
C
      I2 = 0
      DO 880 I = 1, N4
         IF (PARHLD(I)) GO TO 880
         I2 = I2 + 1
         W2(I2) = W2(LEW10-1+I)
  880 CONTINUE
C
C     SET IW(1) TO N2 ,IW(2) TO LEW10 AND W2(4*N2+1),...,W2(4*N2+N4)
C     TO THE 'COMPONENTS' OF PARHLD
C
      IW(1) = N2
      IW(2) = LEW10
      DO 900 I = 1, N4
         IF (PARHLD(I)) THEN
            W2(4*N2+I) = 0.0D0
         ELSE
            W2(4*N2+I) = 1.0D0
         END IF
  900 CONTINUE
C
C     CALL NEW HESSIAN ROUTINE (E04XAF)
C
      I6 = K7 + 1
      DO 920 I = 1, K7
         W2(LEW2-1+I) = 0.0D0
  920 CONTINUE
      IF (K7.GT.0) CALL E04XAF(MSGLVL,K7,EPSA,W2(1),MODE,G13DCS,K7,
     *                         W2(LEW2),F,W2(I6),W2(LEW1),W2(LEW8),
     *                         IWARN,W2(LEW3),IW,W2,IW(LIW2),IDIAG)
C
C     TAKE CARE OF ABNORMAL EXIT FROM FUNCT
C
      IFP = 0
      IF (IDIAG.EQ.2) THEN
         DO 940 I = 1, K7
            IF (IW(LIW2-1+I).GT.0 .AND. IW(LIW2-1+I).LT.4) IFP = 7
            IF (( .NOT. SETQQ) .AND. (IW(LIW2-1+I).EQ.4)) IFP = 7
  940    CONTINUE
      END IF
      IF ((IFP.EQ.7) .AND. (IFAIL.LE.1)) THEN
         WRITE (P01REC,FMT=99999)
         IFAIL = P01ABF(IFAIL,7,SRNAME,2,P01REC)
      END IF
C
C     COPY GRADIENT VECTOR ONTO G
C
      IF ((IDIAG.EQ.-1 .OR. IFP.EQ.7) .AND. (K7.GT.0)) GO TO 980
      I2 = 1
      DO 960 I = 1, N4
         IF (PARHLD(I)) THEN
            G(I) = 0.0D0
         ELSE
            G(I) = -W2(K7+I2)
            I2 = I2 + 1
         END IF
  960 CONTINUE
C
  980 IF (MEAN) THEN
         DO 1000 I = 1, K
            X(K5+I) = X(K5+I) + W2(LW18-1+I)
 1000    CONTINUE
      END IF
C
      IFAILY = 0
C
      IF ((IDIAG.EQ.-1 .OR. IFP.EQ.7) .AND. (K7.GT.0)) THEN
         IFAILY = 6
         DO 1040 I = 1, N4
            DO 1020 J = 1, N4
               DISP(I,J) = 0.0D0
 1020       CONTINUE
            G(I) = 0.0D0
 1040    CONTINUE
         IF (IFAIL.LE.1) THEN
            WRITE (P01REC,FMT=99998)
            IFAIL = P01ABF(IFAIL,IFAILY,SRNAME,2,P01REC)
         END IF
      END IF
C
      IF (IFAILY.EQ.6) GO TO 1140
C
      IF (K7.GT.0) THEN
         IFAILY = 1
         CALL F01AAZ(W2(LEW8),K7,K7,DISP,IDISP,W2(LEW1),IFAILY)
C
         IF (IFAILY.NE.0) THEN
            DO 1080 I = 1, N4
               DO 1060 J = 1, N4
                  DISP(I,J) = 0.0D0
 1060          CONTINUE
 1080       CONTINUE
         END IF
         IF ((IFAILY.NE.0) .AND. (IFAIL.GT.1)) GO TO 1140
         IF ((IFAILY.NE.0) .AND. (IFAIL.LE.1)) THEN
            IFAILY = 8
            WRITE (P01REC,FMT=99993)
     *        '** THE HESSIAN MATRIX IS NOT POSITIVE-DEFINITE',
     *        ' AT THE SOLUTION POINT'
            IFAIL = P01ABF(IFAIL,IFAILY,SRNAME,1,P01REC)
            GO TO 1140
         END IF
      END IF
C
C     PUT ZERO ROWS AND COLUMNS INTO DISP
C
      IF (IBOUND.EQ.0) THEN
         CALL G13DCZ(W2(LEW8),N4,PARHLD,DISP,IDISP)
         DO 1120 I = 1, N4
            DO 1100 J = 1, N4
               DISP(I,J) = W2(LEW8-1+(J-1)*N4+I)
 1100       CONTINUE
 1120    CONTINUE
      END IF
C
 1140 IF (IFAILY.NE.0) GO TO 1300
C
C     IF HESSIAN MATRIX WAS SUCCESSFULLY INVERTED THEN CORRECT
C     COVARIANCE MATRIX TO CORRELATION - STANDARD ERROR MATRIX
C
      DO 1240 I = 1, N4
         IF ( .NOT. PARHLD(I)) THEN
            EPSA = DISP(I,I)
            IF (EPSA.LE.0.0D0) THEN
               IFAILY = 8
               DO 1180 ITT = 1, N4
                  DO 1160 J = 1, N4
                     DISP(ITT,J) = 0.0D0
 1160             CONTINUE
 1180          CONTINUE
               IF (IFAIL.LE.1) IFAIL = P01ABF(IFAIL,IFAILY,SRNAME,0,
     *                                 P01REC)
               GO TO 1300
            END IF
            DISP(I,I) = SQRT(DISP(I,I))
            IF (I.EQ.1) GO TO 1240
            DO 1200 J = 1, I - 1
               IF ( .NOT. PARHLD(J)) THEN
                  DISP(I,J) = DISP(I,J)/SQRT(EPSA*DISP(J,J))
               ELSE
                  DISP(I,J) = 0.0D0
               END IF
 1200       CONTINUE
         ELSE
            DO 1220 J = 1, I
               DISP(I,J) = 0.0D0
 1220       CONTINUE
         END IF
 1240 CONTINUE
C
C     ASSIGN UPPER TRIANGLE OF DISP
C
      DO 1280 I = 1, N4
         DO 1260 J = 1, I
            DISP(J,I) = DISP(I,J)
 1260    CONTINUE
 1280 CONTINUE
C
C     IF ISHOW = 1 OR 2 DISPLAY PARAMETER ESTIMATES AND THEIR
C     STANDARD ERRORS (PLUS RESIDUALS)
C
 1300 IF (IFAIL.LE.1) IFAIL = 0
      IF (ISHOW.GT.0) CALL G13DCT(LP,IFAIL,LOGL,P,Q,K,X,N4,DISP,IDISP,
     *                            MEAN,QQ,IK,ISHOW,N,V)
C
      RETURN
C
99999 FORMAT (' ** THE ELEMENTS OF THE 2 ND DERIVATIVE MATRIX AND THE',
     *       /'    GRADIENT VECTOR COULD NOT BE EVALUATED ACCURATELY')
99998 FORMAT (' ** THE SOLUTION IS TOO CLOSE TO THE ADMISSIBILITY BOUN',
     *       'DARY FOR STANDARD ERRORS',/'    AND CORRELATIONS OF PARA',
     *       'METER ESTIMATES TO BE CALCULATED')
99997 FORMAT ('  ** ON ENTRY, THE VALUE OF CGETOL ( ',1P,D12.5,' ) IS ',
     *       'TOO SMALL',/' AND THE VALUE ',1P,D12.5,' HAS BEEN USED I',
     *       'NSTEAD')
99996 FORMAT (' ** ON ENTRY, LWORK ( ',I16,' ) MUST BE AT LEAST ',I16)
99995 FORMAT (' ** ON ENTRY, LIW ( ',I16,' ) MUST BE AT LEAST ',I16)
99994 FORMAT (1X,A)
99993 FORMAT (1X,A,A)
99992 FORMAT (1X,A,I16,A)
99991 FORMAT (' ** ON ENTRY, ONE OR MORE OF THE FOLLOWING PARAMETER VA',
     *       'LUES IS ILLEGAL',/' K      = ',I16,'   N     = ',I16,'  ',
     *       ' IP  = ',I16,/' IQ     = ',I16,'   NPAR  = ',I16,'   IK ',
     *       ' = ',I16,/' MAXCAL = ',I16,'   ISHOW = ',I16,'   ICM = ',
     *       I16,/' MEAN   = ',L16)
      END
