      SUBROUTINE E04FCY(M,N,LH,NS,EPSMCH,LSFJSA,LSF,IGRADE,X,FVEC,FF,
     *                  P1ZERO,V,IV,P,PHESL,PHESD,RHS,IFLAG,IW,LIW,W,LW)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C
C     **************************************************************
C
C     THIS SUBROUTINE COMPUTES FINITE-DIFFERENCE APPROXIMATIONS TO
C     THE MATRIX VTHV, WHERE H IS THE SECOND TERM OF THE HESSIAN
C     MATRIX OF THE SUM OF SQUARES, AND THE VECTOR VTHP.
C     THE WORKSPACE ARRAY W IS OF LENGTH AT LEAST 2*N + M.
C
C     PHILIP E. GILL, ENID M. R. LONG, WALTER MURRAY,
C     SUSAN M. PICKEN AND BRIAN T. HINDE
C     D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND
C
C     **************************************************************
C
C     Modified to call BLAS.
C     Peter Mayes, NAG Central Office, October 1987.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPSMCH, FF
      INTEGER           IFLAG, IGRADE, IV, LH, LIW, LW, M, N, NS
      LOGICAL           P1ZERO
C     .. Array Arguments ..
      DOUBLE PRECISION  FVEC(M), P(N), PHESD(NS), PHESL(LH), RHS(NS),
     *                  V(IV,N), W(LW), X(N)
      INTEGER           IW(LIW)
C     .. Subroutine Arguments ..
      EXTERNAL          LSF, LSFJSA
C     .. Local Scalars ..
      DOUBLE PRECISION  FF1HV, FF2HV, H, H2, H2INV, HP, PNORM, RHSJ,
     *                  STEP, X1HVK
      INTEGER           DUMMY, I, IR, J, JR, JUMP, K, K1, K2, L, LCOL,
     *                  LD, LFPLUS, LPHESL, LROW, LX1HV, LX2HV
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2
      EXTERNAL          DDOT, DNRM2
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
C
C     ALLOCATE THE WORKSPACE AND COMPUTE THE FINITE-DIFFERENCE
C     INTERVALS.
C
      IF (NS.LE.0) RETURN
      DUMMY = 1
      LD = M
      LX1HV = 1
      LX2HV = LX1HV + N
      LFPLUS = LX2HV + N
      H2 = SQRT(EPSMCH)
      H = SQRT(H2)
      H2INV = 1.0D+0/H2
      IF (P1ZERO) GO TO 20
      PNORM = DNRM2(N,P,1)
      HP = H/PNORM
C
C     INITIALIZE THE NON-DIAGONAL ELEMENTS OF VTHV.
C
   20 IF (NS.EQ.1) GO TO 60
      LPHESL = NS*(NS-1)/2
      DO 40 L = 1, LPHESL
         PHESL(L) = FF
   40 CONTINUE
   60 DO 220 J = 1, NS
C
C        COMPUTE THE FUNCTION VALUES AT INTERVALS OF H AND 2*H ALONG
C        THE JTH COLUMN OF V.
C
         JR = IGRADE + J
         K1 = LX1HV
         K2 = LX2HV
         DO 80 K = 1, N
            STEP = H*V(K,JR)
            X1HVK = X(K) + STEP
            W(K1) = X1HVK
            W(K2) = X1HVK + STEP
            K1 = K1 + 1
            K2 = K2 + 1
   80    CONTINUE
         IFLAG = 0
         CALL LSFJSA(IFLAG,M,N,LSF,W(LX1HV),W(LFPLUS),W(DUMMY),LD,IW,
     *               LIW,W,LW)
         IF (IFLAG.LT.0) RETURN
         FF1HV = DDOT(M,FVEC,1,W(LFPLUS),1)
         IFLAG = 0
         CALL LSFJSA(IFLAG,M,N,LSF,W(LX2HV),W(LFPLUS),W(DUMMY),LD,IW,
     *               LIW,W,LW)
         IF (IFLAG.LT.0) RETURN
         FF2HV = DDOT(M,FVEC,1,W(LFPLUS),1)
C
C        COMPUTE THE JTH DIAGONAL ELEMENT OF VTHV AND INITIALIZE
C        THE JTH ELEMENT OF VTHP.
C
         PHESD(J) = (FF2HV-FF1HV-FF1HV+FF)*H2INV
         RHSJ = 0.0D+0
         IF (P1ZERO) GO TO 120
         K1 = LX1HV
         K2 = LX2HV
         DO 100 K = 1, N
            W(K2) = W(K1) + HP*P(K)
            K1 = K1 + 1
            K2 = K2 + 1
  100    CONTINUE
         IFLAG = 0
         CALL LSFJSA(IFLAG,M,N,LSF,W(LX2HV),W(LFPLUS),W(DUMMY),LD,IW,
     *               LIW,W,LW)
         IF (IFLAG.LT.0) RETURN
         FF2HV = DDOT(M,FVEC,1,W(LFPLUS),1)
         RHSJ = FF2HV - FF1HV + FF
  120    RHS(J) = RHSJ
C
C        UPDATE THE NON-DIAGONAL ELEMENTS OF THE JTH ROW OF VTHV.
C
         IF (J.EQ.1) GO TO 180
         LROW = J - 1
         L = LROW*(LROW-1)/2
         DO 160 I = 1, LROW
            IR = IGRADE + I
            K1 = LX1HV
            K2 = LX2HV
            DO 140 K = 1, N
               W(K2) = W(K1) + H*V(K,IR)
               K1 = K1 + 1
               K2 = K2 + 1
  140       CONTINUE
            IFLAG = 0
            CALL LSFJSA(IFLAG,M,N,LSF,W(LX2HV),W(LFPLUS),W(DUMMY),LD,IW,
     *                  LIW,W,LW)
            IF (IFLAG.LT.0) RETURN
            FF2HV = DDOT(M,FVEC,1,W(LFPLUS),1)
            L = L + 1
            PHESL(L) = PHESL(L) - FF1HV + FF2HV
  160    CONTINUE
C
C        UPDATE THE NON-DIAGONAL ELEMENTS OF THE JTH COLUMN OF VTHV.
C
  180    IF (J.EQ.NS) GO TO 220
         JUMP = J
         LCOL = NS - J
         L = J*(J+1)/2
         DO 200 I = 1, LCOL
            PHESL(L) = PHESL(L) - FF1HV
            L = L + JUMP
            JUMP = JUMP + 1
  200    CONTINUE
  220 CONTINUE
C
C     COMPLETE THE COMPUTATION OF THE NON-DIAGONAL ELEMENTS OF VTHV.
C
      IF (NS.EQ.1) GO TO 260
      DO 240 L = 1, LPHESL
         PHESL(L) = PHESL(L)*H2INV
  240 CONTINUE
C
C     COMPLETE THE COMPUTATION OF THE VECTOR VTHP.
C
  260 IF (P1ZERO) RETURN
      K1 = LX1HV
      DO 280 K = 1, N
         W(K1) = X(K) + HP*P(K)
         K1 = K1 + 1
  280 CONTINUE
      IFLAG = 0
      CALL LSFJSA(IFLAG,M,N,LSF,W(LX1HV),W(LFPLUS),W(DUMMY),LD,IW,LIW,W,
     *            LW)
      IF (IFLAG.LT.0) RETURN
      FF1HV = DDOT(M,FVEC,1,W(LFPLUS),1)
      DO 300 J = 1, NS
         RHS(J) = (RHS(J)-FF1HV)*H2INV*PNORM
  300 CONTINUE
      RETURN
C
C     END OF E04FCY   (FDVTHV)
C
      END
