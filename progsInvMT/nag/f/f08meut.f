      SUBROUTINE F08MEU(N,Q,E,QQ,EE,SUP,SIGMA,KEND,OFF,IPHASE,ICONV,EPS,
     *                  TOL2,SMALL2)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Purpose
C     =======
C
C     DLASQ3 is the workhorse of the whole bidiagonal SVD
C     algorithm. This can be described as the differential qd
C     with shifts.
C
C     Arguments
C     =========
C
C  N       (input/output) INTEGER
C          On entry, N specifies the number of rows and columns
C          in the matrix. N must be at least 3.
C          On exit N is non-negative and less than the input value.
C
C  Q       (input/output) DOUBLE PRECISION array, dimension (N)
C          Q array in ping (see IPHASE below)
C
C  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
C          E array in ping (see IPHASE below)
C
C  QQ      (input/output) DOUBLE PRECISION array, dimension (N)
C          Q array in pong (see IPHASE below)
C
C  EE      (input/output) DOUBLE PRECISION array, dimension (N-1)
C          E array in pong (see IPHASE below)
C
C  SUP     (input/output) DOUBLE PRECISION
C          Upper bound for the smallest eigenvalue
C
C  SIGMA   (input/output) DOUBLE PRECISION
C          Accumulated shift for the present submatrix
C
C  KEND    (input/output) INTEGER
C          Index where minimum D(i) occurs in recurrence for
C          splitting criterion
C
C  OFF     (input/output) INTEGER
C          Offset for arrays
C
C  IPHASE  (input/output) INTEGER
C          If IPHASE = 1 (ping) then data is in Q and E arrays
C          If IPHASE = 2 (pong) then data is in QQ and EE arrays
C
C  ICONV   (input) INTEGER
C          If ICONV = 0 a bottom part of a matrix (with a split)
C          If ICONV =-3 a top part of a matrix (with a split)
C
C  EPS     (input) DOUBLE PRECISION
C          Machine epsilon
C
C  TOL2    (input) DOUBLE PRECISION
C          Square of the relative tolerance TOL as defined in DLASQ1
C
C  SMALL2  (input) DOUBLE PRECISION
C          A threshold value as defined in DLASQ1
C
C  -- LAPACK routine (version 2.0) (adapted NAG library) --
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C     September 30, 1994
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
      INTEGER           NPP
      PARAMETER         (NPP=32)
      INTEGER           IPP
      PARAMETER         (IPP=5)
      DOUBLE PRECISION  HALF, FOUR
      PARAMETER         (HALF=0.5D+0,FOUR=4.0D+0)
      INTEGER           IFLMAX
      PARAMETER         (IFLMAX=2)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS, SIGMA, SMALL2, SUP, TOL2
      INTEGER           ICONV, IPHASE, KEND, N, OFF
C     .. Array Arguments ..
      DOUBLE PRECISION  E(*), EE(*), Q(*), QQ(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  D, DM, QEMAX, T1, TAU, TOLX, TOLY, TOLZ, XX, YY
      INTEGER           I, IC, ICNT, IFL, IP, ISP, K1END, K2END, KE, KS,
     *                  MAXIT, N1, N2
      LOGICAL           LDEF, LSPLIT
C     .. External Subroutines ..
      EXTERNAL          DCOPY, F08MET
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Executable Statements ..
      ICNT = 0
      TAU = ZERO
      DM = SUP
      TOLX = SIGMA*TOL2
      TOLZ = MAX(SMALL2,SIGMA)*TOL2
C
C     Set maximum number of iterations
C
      MAXIT = 100*N
C
C     Flipping
C
      IC = 2
      IF (N.GT.3) THEN
         IF (IPHASE.EQ.1) THEN
            DO 20 I = 1, N - 2
               IF (Q(I).GT.Q(I+1)) IC = IC + 1
               IF (E(I).GT.E(I+1)) IC = IC + 1
   20       CONTINUE
            IF (Q(N-1).GT.Q(N)) IC = IC + 1
            IF (IC.LT.N) THEN
               CALL DCOPY(N,Q,1,QQ,-1)
               CALL DCOPY(N-1,E,1,EE,-1)
               IF (KEND.NE.0) KEND = N - KEND + 1
               IPHASE = 2
            END IF
         ELSE
            DO 40 I = 1, N - 2
               IF (QQ(I).GT.QQ(I+1)) IC = IC + 1
               IF (EE(I).GT.EE(I+1)) IC = IC + 1
   40       CONTINUE
            IF (QQ(N-1).GT.QQ(N)) IC = IC + 1
            IF (IC.LT.N) THEN
               CALL DCOPY(N,QQ,1,Q,-1)
               CALL DCOPY(N-1,EE,1,E,-1)
               IF (KEND.NE.0) KEND = N - KEND + 1
               IPHASE = 1
            END IF
         END IF
      END IF
      IF (ICONV.EQ.-3) THEN
         IF (IPHASE.EQ.1) THEN
            GO TO 360
         ELSE
            GO TO 160
         END IF
      END IF
      IF (IPHASE.EQ.2) GO TO 260
C
C     The ping section of the code
C
   60 CONTINUE
      IFL = 0
C
C     Compute the shift
C
      IF (KEND.EQ.0 .OR. SUP.EQ.ZERO) THEN
         TAU = ZERO
      ELSE IF (ICNT.GT.0 .AND. DM.LE.TOLZ) THEN
         TAU = ZERO
      ELSE
         IP = MAX(IPP,N/NPP)
         N2 = 2*IP + 1
         IF (N2.GE.N) THEN
            N1 = 1
            N2 = N
         ELSE IF (KEND+IP.GT.N) THEN
            N1 = N - 2*IP
         ELSE IF (KEND-IP.LT.1) THEN
            N1 = 1
         ELSE
            N1 = KEND - IP
         END IF
         CALL F08MET(N2,Q(N1),E(N1),TAU,SUP)
      END IF
   80 CONTINUE
      ICNT = ICNT + 1
      IF (ICNT.GT.MAXIT) THEN
         SUP = -ONE
         RETURN
      END IF
      IF (TAU.EQ.ZERO) THEN
C
C     dqd algorithm
C
         D = Q(1)
         DM = D
         KE = 0
         DO 100 I = 1, N - 3
            QQ(I) = D + E(I)
            D = (D/QQ(I))*Q(I+1)
            IF (DM.GT.D) THEN
               DM = D
               KE = I
            END IF
  100    CONTINUE
         KE = KE + 1
C
C     Penultimate dqd step (in ping)
C
         K2END = KE
         QQ(N-2) = D + E(N-2)
         D = (D/QQ(N-2))*Q(N-1)
         IF (DM.GT.D) THEN
            DM = D
            KE = N - 1
         END IF
C
C     Final dqd step (in ping)
C
         K1END = KE
         QQ(N-1) = D + E(N-1)
         D = (D/QQ(N-1))*Q(N)
         IF (DM.GT.D) THEN
            DM = D
            KE = N
         END IF
         QQ(N) = D
      ELSE
C
C     The dqds algorithm (in ping)
C
         D = Q(1) - TAU
         DM = D
         KE = 0
         IF (D.LT.ZERO) GO TO 240
         DO 120 I = 1, N - 3
            QQ(I) = D + E(I)
            D = (D/QQ(I))*Q(I+1) - TAU
            IF (DM.GT.D) THEN
               DM = D
               KE = I
               IF (D.LT.ZERO) GO TO 240
            END IF
  120    CONTINUE
         KE = KE + 1
C
C     Penultimate dqds step (in ping)
C
         K2END = KE
         QQ(N-2) = D + E(N-2)
         D = (D/QQ(N-2))*Q(N-1) - TAU
         IF (DM.GT.D) THEN
            DM = D
            KE = N - 1
            IF (D.LT.ZERO) GO TO 240
         END IF
C
C     Final dqds step (in ping)
C
         K1END = KE
         QQ(N-1) = D + E(N-1)
         D = (D/QQ(N-1))*Q(N) - TAU
         IF (DM.GT.D) THEN
            DM = D
            KE = N
         END IF
         QQ(N) = D
      END IF
C
C        Convergence when QQ(N) is small (in ping)
C
      IF (ABS(QQ(N)).LE.SIGMA*TOL2) THEN
         QQ(N) = ZERO
         DM = ZERO
         KE = N
      END IF
      IF (QQ(N).LT.ZERO) GO TO 240
C
C     Non-negative qd array: Update the e's
C
      DO 140 I = 1, N - 1
         EE(I) = (E(I)/QQ(I))*Q(I+1)
  140 CONTINUE
C
C     Updating sigma and iphase in ping
C
      SIGMA = SIGMA + TAU
      IPHASE = 2
  160 CONTINUE
      TOLX = SIGMA*TOL2
      TOLY = SIGMA*EPS
      TOLZ = MAX(SIGMA,SMALL2)*TOL2
C
C     Checking for deflation and convergence (in ping)
C
  180 CONTINUE
      IF (N.LE.2) RETURN
C
C        Deflation: bottom 1x1 (in ping)
C
      LDEF = .FALSE.
      IF (EE(N-1).LE.TOLZ) THEN
         LDEF = .TRUE.
      ELSE IF (SIGMA.GT.ZERO) THEN
         IF (EE(N-1).LE.EPS*(SIGMA+QQ(N))) THEN
            IF (EE(N-1)*(QQ(N)/(QQ(N)+SIGMA)).LE.TOL2*(QQ(N)+SIGMA))
     *          THEN
               LDEF = .TRUE.
            END IF
         END IF
      ELSE
         IF (EE(N-1).LE.QQ(N)*TOL2) THEN
            LDEF = .TRUE.
         END IF
      END IF
      IF (LDEF) THEN
         Q(N) = QQ(N) + SIGMA
         N = N - 1
         ICONV = ICONV + 1
         GO TO 180
      END IF
C
C        Deflation: bottom 2x2 (in ping)
C
      LDEF = .FALSE.
      IF (EE(N-2).LE.TOLZ) THEN
         LDEF = .TRUE.
      ELSE IF (SIGMA.GT.ZERO) THEN
         T1 = SIGMA + EE(N-1)*(SIGMA/(SIGMA+QQ(N)))
         IF (EE(N-2)*(T1/(QQ(N-1)+T1)).LE.TOLY) THEN
            IF (EE(N-2)*(QQ(N-1)/(QQ(N-1)+T1)).LE.TOLX) THEN
               LDEF = .TRUE.
            END IF
         END IF
      ELSE
         IF (EE(N-2).LE.(QQ(N)/(QQ(N)+EE(N-1)+QQ(N-1)))*QQ(N-1)*TOL2)
     *       THEN
            LDEF = .TRUE.
         END IF
      END IF
      IF (LDEF) THEN
         QEMAX = MAX(QQ(N),QQ(N-1),EE(N-1))
         IF (QEMAX.NE.ZERO) THEN
            IF (QEMAX.EQ.QQ(N-1)) THEN
               XX = HALF*(QQ(N)+QQ(N-1)+EE(N-1)+QEMAX*SQRT(((QQ(N)
     *              -QQ(N-1)+EE(N-1))/QEMAX)**2+FOUR*EE(N-1)/QEMAX))
            ELSE IF (QEMAX.EQ.QQ(N)) THEN
               XX = HALF*(QQ(N)+QQ(N-1)+EE(N-1)+QEMAX*SQRT(((QQ(N-1)
     *              -QQ(N)+EE(N-1))/QEMAX)**2+FOUR*EE(N-1)/QEMAX))
            ELSE
               XX = HALF*(QQ(N)+QQ(N-1)+EE(N-1)+QEMAX*SQRT(((QQ(N)
     *              -QQ(N-1)+EE(N-1))/QEMAX)**2+FOUR*QQ(N-1)/QEMAX))
            END IF
            YY = (MAX(QQ(N),QQ(N-1))/XX)*MIN(QQ(N),QQ(N-1))
         ELSE
            XX = ZERO
            YY = ZERO
         END IF
         Q(N-1) = SIGMA + XX
         Q(N) = YY + SIGMA
         N = N - 2
         ICONV = ICONV + 2
         GO TO 180
      END IF
C
C     Updating bounds before going to pong
C
      IF (ICONV.EQ.0) THEN
         KEND = KE
         SUP = MIN(DM,SUP-TAU)
      ELSE IF (ICONV.GT.0) THEN
         SUP = MIN(QQ(N),QQ(N-1),QQ(N-2),QQ(1),QQ(2),QQ(3))
         IF (ICONV.EQ.1) THEN
            KEND = K1END
         ELSE IF (ICONV.EQ.2) THEN
            KEND = K2END
         ELSE
            KEND = N
         END IF
         ICNT = 0
         MAXIT = 100*N
      END IF
C
C     Checking for splitting in ping
C
      LSPLIT = .FALSE.
      DO 200 KS = N - 3, 3, -1
         IF (EE(KS).LE.TOLY) THEN
            IF (EE(KS)*(MIN(QQ(KS+1),QQ(KS))/(MIN(QQ(KS+1),QQ(KS))
     *          +SIGMA)).LE.TOLX) THEN
               LSPLIT = .TRUE.
               GO TO 220
            END IF
         END IF
  200 CONTINUE
C
      KS = 2
      IF (EE(2).LE.TOLZ) THEN
         LSPLIT = .TRUE.
      ELSE IF (SIGMA.GT.ZERO) THEN
         T1 = SIGMA + EE(1)*(SIGMA/(SIGMA+QQ(1)))
         IF (EE(2)*(T1/(QQ(1)+T1)).LE.TOLY) THEN
            IF (EE(2)*(QQ(1)/(QQ(1)+T1)).LE.TOLX) THEN
               LSPLIT = .TRUE.
            END IF
         END IF
      ELSE
         IF (EE(2).LE.(QQ(1)/(QQ(1)+EE(1)+QQ(2)))*QQ(2)*TOL2) THEN
            LSPLIT = .TRUE.
         END IF
      END IF
      IF (LSPLIT) GO TO 220
C
      KS = 1
      IF (EE(1).LE.TOLZ) THEN
         LSPLIT = .TRUE.
      ELSE IF (SIGMA.GT.ZERO) THEN
         IF (EE(1).LE.EPS*(SIGMA+QQ(1))) THEN
            IF (EE(1)*(QQ(1)/(QQ(1)+SIGMA)).LE.TOL2*(QQ(1)+SIGMA)) THEN
               LSPLIT = .TRUE.
            END IF
         END IF
      ELSE
         IF (EE(1).LE.QQ(1)*TOL2) THEN
            LSPLIT = .TRUE.
         END IF
      END IF
C
  220 CONTINUE
      IF (LSPLIT) THEN
         SUP = MIN(QQ(N),QQ(N-1),QQ(N-2))
         ISP = -(OFF+1)
         OFF = OFF + KS
         N = N - KS
         KEND = MAX(1,KEND-KS)
         E(KS) = SIGMA
         EE(KS) = ISP
         ICONV = 0
         RETURN
      END IF
C
C     Coincidence
C
      IF (TAU.EQ.ZERO .AND. DM.LE.TOLZ .AND. KEND.NE.N .AND. ICONV.EQ.
     *    0 .AND. ICNT.GT.0) THEN
         CALL DCOPY(N-KE,E(KE),1,QQ(KE),1)
         QQ(N) = ZERO
         CALL DCOPY(N-KE,Q(KE+1),1,EE(KE),1)
         SUP = ZERO
      END IF
      ICONV = 0
      GO TO 260
C
C     A new shift when the previous failed (in ping)
C
  240 CONTINUE
      IFL = IFL + 1
      SUP = TAU
C
C     SUP is small or
C     Too many bad shifts (ping)
C
      IF (SUP.LE.TOLZ .OR. IFL.GE.IFLMAX) THEN
         TAU = ZERO
         GO TO 80
C
C     The asymptotic shift (in ping)
C
      ELSE
         TAU = MAX(TAU+D,ZERO)
         IF (TAU.LE.TOLZ) TAU = ZERO
         GO TO 80
      END IF
C
C     the pong section of the code
C
  260 CONTINUE
      IFL = 0
C
C     Compute the shift (in pong)
C
      IF (KEND.EQ.0 .AND. SUP.EQ.ZERO) THEN
         TAU = ZERO
      ELSE IF (ICNT.GT.0 .AND. DM.LE.TOLZ) THEN
         TAU = ZERO
      ELSE
         IP = MAX(IPP,N/NPP)
         N2 = 2*IP + 1
         IF (N2.GE.N) THEN
            N1 = 1
            N2 = N
         ELSE IF (KEND+IP.GT.N) THEN
            N1 = N - 2*IP
         ELSE IF (KEND-IP.LT.1) THEN
            N1 = 1
         ELSE
            N1 = KEND - IP
         END IF
         CALL F08MET(N2,QQ(N1),EE(N1),TAU,SUP)
      END IF
  280 CONTINUE
      ICNT = ICNT + 1
      IF (ICNT.GT.MAXIT) THEN
         SUP = -SUP
         RETURN
      END IF
      IF (TAU.EQ.ZERO) THEN
C
C     The dqd algorithm (in pong)
C
         D = QQ(1)
         DM = D
         KE = 0
         DO 300 I = 1, N - 3
            Q(I) = D + EE(I)
            D = (D/Q(I))*QQ(I+1)
            IF (DM.GT.D) THEN
               DM = D
               KE = I
            END IF
  300    CONTINUE
         KE = KE + 1
C
C     Penultimate dqd step (in pong)
C
         K2END = KE
         Q(N-2) = D + EE(N-2)
         D = (D/Q(N-2))*QQ(N-1)
         IF (DM.GT.D) THEN
            DM = D
            KE = N - 1
         END IF
C
C     Final dqd step (in pong)
C
         K1END = KE
         Q(N-1) = D + EE(N-1)
         D = (D/Q(N-1))*QQ(N)
         IF (DM.GT.D) THEN
            DM = D
            KE = N
         END IF
         Q(N) = D
      ELSE
C
C     The dqds algorithm (in pong)
C
         D = QQ(1) - TAU
         DM = D
         KE = 0
         IF (D.LT.ZERO) GO TO 440
         DO 320 I = 1, N - 3
            Q(I) = D + EE(I)
            D = (D/Q(I))*QQ(I+1) - TAU
            IF (DM.GT.D) THEN
               DM = D
               KE = I
               IF (D.LT.ZERO) GO TO 440
            END IF
  320    CONTINUE
         KE = KE + 1
C
C     Penultimate dqds step (in pong)
C
         K2END = KE
         Q(N-2) = D + EE(N-2)
         D = (D/Q(N-2))*QQ(N-1) - TAU
         IF (DM.GT.D) THEN
            DM = D
            KE = N - 1
            IF (D.LT.ZERO) GO TO 440
         END IF
C
C     Final dqds step (in pong)
C
         K1END = KE
         Q(N-1) = D + EE(N-1)
         D = (D/Q(N-1))*QQ(N) - TAU
         IF (DM.GT.D) THEN
            DM = D
            KE = N
         END IF
         Q(N) = D
      END IF
C
C        Convergence when Q(N) is small (in pong)
C
      IF (ABS(Q(N)).LE.SIGMA*TOL2) THEN
         Q(N) = ZERO
         DM = ZERO
         KE = N
      END IF
      IF (Q(N).LT.ZERO) GO TO 440
C
C     Non-negative qd array: Update the e's
C
      DO 340 I = 1, N - 1
         E(I) = (EE(I)/Q(I))*QQ(I+1)
  340 CONTINUE
C
C     Updating sigma and iphase in pong
C
      SIGMA = SIGMA + TAU
  360 CONTINUE
      IPHASE = 1
      TOLX = SIGMA*TOL2
      TOLY = SIGMA*EPS
C
C     Checking for deflation and convergence (in pong)
C
  380 CONTINUE
      IF (N.LE.2) RETURN
C
C        Deflation: bottom 1x1 (in pong)
C
      LDEF = .FALSE.
      IF (E(N-1).LE.TOLZ) THEN
         LDEF = .TRUE.
      ELSE IF (SIGMA.GT.ZERO) THEN
         IF (E(N-1).LE.EPS*(SIGMA+Q(N))) THEN
            IF (E(N-1)*(Q(N)/(Q(N)+SIGMA)).LE.TOL2*(Q(N)+SIGMA)) THEN
               LDEF = .TRUE.
            END IF
         END IF
      ELSE
         IF (E(N-1).LE.Q(N)*TOL2) THEN
            LDEF = .TRUE.
         END IF
      END IF
      IF (LDEF) THEN
         Q(N) = Q(N) + SIGMA
         N = N - 1
         ICONV = ICONV + 1
         GO TO 380
      END IF
C
C        Deflation: bottom 2x2 (in pong)
C
      LDEF = .FALSE.
      IF (E(N-2).LE.TOLZ) THEN
         LDEF = .TRUE.
      ELSE IF (SIGMA.GT.ZERO) THEN
         T1 = SIGMA + E(N-1)*(SIGMA/(SIGMA+Q(N)))
         IF (E(N-2)*(T1/(Q(N-1)+T1)).LE.TOLY) THEN
            IF (E(N-2)*(Q(N-1)/(Q(N-1)+T1)).LE.TOLX) THEN
               LDEF = .TRUE.
            END IF
         END IF
      ELSE
         IF (E(N-2).LE.(Q(N)/(Q(N)+EE(N-1)+Q(N-1))*Q(N-1))*TOL2) THEN
            LDEF = .TRUE.
         END IF
      END IF
      IF (LDEF) THEN
         QEMAX = MAX(Q(N),Q(N-1),E(N-1))
         IF (QEMAX.NE.ZERO) THEN
            IF (QEMAX.EQ.Q(N-1)) THEN
               XX = HALF*(Q(N)+Q(N-1)+E(N-1)+QEMAX*SQRT(((Q(N)-Q(N-1)
     *              +E(N-1))/QEMAX)**2+FOUR*E(N-1)/QEMAX))
            ELSE IF (QEMAX.EQ.Q(N)) THEN
               XX = HALF*(Q(N)+Q(N-1)+E(N-1)+QEMAX*SQRT(((Q(N-1)-Q(N)
     *              +E(N-1))/QEMAX)**2+FOUR*E(N-1)/QEMAX))
            ELSE
               XX = HALF*(Q(N)+Q(N-1)+E(N-1)+QEMAX*SQRT(((Q(N)-Q(N-1)
     *              +E(N-1))/QEMAX)**2+FOUR*Q(N-1)/QEMAX))
            END IF
            YY = (MAX(Q(N),Q(N-1))/XX)*MIN(Q(N),Q(N-1))
         ELSE
            XX = ZERO
            YY = ZERO
         END IF
         Q(N-1) = SIGMA + XX
         Q(N) = YY + SIGMA
         N = N - 2
         ICONV = ICONV + 2
         GO TO 380
      END IF
C
C     Updating bounds before going to pong
C
      IF (ICONV.EQ.0) THEN
         KEND = KE
         SUP = MIN(DM,SUP-TAU)
      ELSE IF (ICONV.GT.0) THEN
         SUP = MIN(Q(N),Q(N-1),Q(N-2),Q(1),Q(2),Q(3))
         IF (ICONV.EQ.1) THEN
            KEND = K1END
         ELSE IF (ICONV.EQ.2) THEN
            KEND = K2END
         ELSE
            KEND = N
         END IF
         ICNT = 0
         MAXIT = 100*N
      END IF
C
C     Checking for splitting in pong
C
      LSPLIT = .FALSE.
      DO 400 KS = N - 3, 3, -1
         IF (E(KS).LE.TOLY) THEN
            IF (E(KS)*(MIN(Q(KS+1),Q(KS))/(MIN(Q(KS+1),Q(KS))+SIGMA))
     *          .LE.TOLX) THEN
               LSPLIT = .TRUE.
               GO TO 420
            END IF
         END IF
  400 CONTINUE
C
      KS = 2
      IF (E(2).LE.TOLZ) THEN
         LSPLIT = .TRUE.
      ELSE IF (SIGMA.GT.ZERO) THEN
         T1 = SIGMA + E(1)*(SIGMA/(SIGMA+Q(1)))
         IF (E(2)*(T1/(Q(1)+T1)).LE.TOLY) THEN
            IF (E(2)*(Q(1)/(Q(1)+T1)).LE.TOLX) THEN
               LSPLIT = .TRUE.
            END IF
         END IF
      ELSE
         IF (E(2).LE.(Q(1)/(Q(1)+E(1)+Q(2)))*Q(2)*TOL2) THEN
            LSPLIT = .TRUE.
         END IF
      END IF
      IF (LSPLIT) GO TO 420
C
      KS = 1
      IF (E(1).LE.TOLZ) THEN
         LSPLIT = .TRUE.
      ELSE IF (SIGMA.GT.ZERO) THEN
         IF (E(1).LE.EPS*(SIGMA+Q(1))) THEN
            IF (E(1)*(Q(1)/(Q(1)+SIGMA)).LE.TOL2*(Q(1)+SIGMA)) THEN
               LSPLIT = .TRUE.
            END IF
         END IF
      ELSE
         IF (E(1).LE.Q(1)*TOL2) THEN
            LSPLIT = .TRUE.
         END IF
      END IF
C
  420 CONTINUE
      IF (LSPLIT) THEN
         SUP = MIN(Q(N),Q(N-1),Q(N-2))
         ISP = OFF + 1
         OFF = OFF + KS
         KEND = MAX(1,KEND-KS)
         N = N - KS
         E(KS) = SIGMA
         EE(KS) = ISP
         ICONV = 0
         RETURN
      END IF
C
C     Coincidence
C
      IF (TAU.EQ.ZERO .AND. DM.LE.TOLZ .AND. KEND.NE.N .AND. ICONV.EQ.
     *    0 .AND. ICNT.GT.0) THEN
         CALL DCOPY(N-KE,EE(KE),1,Q(KE),1)
         Q(N) = ZERO
         CALL DCOPY(N-KE,QQ(KE+1),1,E(KE),1)
         SUP = ZERO
      END IF
      ICONV = 0
      GO TO 60
C
C     Computation of a new shift when the previous failed (in pong)
C
  440 CONTINUE
      IFL = IFL + 1
      SUP = TAU
C
C     SUP is small or
C     Too many bad shifts (in pong)
C
      IF (SUP.LE.TOLZ .OR. IFL.GE.IFLMAX) THEN
         TAU = ZERO
         GO TO 280
C
C     The asymptotic shift (in pong)
C
      ELSE
         TAU = MAX(TAU+D,ZERO)
         IF (TAU.LE.TOLZ) TAU = ZERO
         GO TO 280
      END IF
C
C     End of F08MEU (DLASQ3)
C
      END
