      SUBROUTINE E04HEW(M,N,LH,NS,EPSMCH,LSFJSA,LSF,IGRADE,X,FVEC,FJAC,
     *                  LJ,V,IV,P,PHESL,PHESD,RHS,IFLAG,IW,LIW,W,LW)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C
C     **************************************************************
C
C     THIS SUBROUTINE COMPUTES THE FINITE-DIFFERENCE APPROXIMATION
C     TO THE PROJECTION OF THE SECOND TERM OF THE HESSIAN MATRIX OF
C     THE SUM OF SQUARES, H. IT ALSO FORMS VTHP, WHERE V AND P ARE
C     GIVEN ARRAYS.
C     THE WORKSPACE ARRAY W IS OF LENGTH AT LEAST 2*N + M*N + M.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN AND
C     NICHOLAS I. M. GOULD.
C     D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND.
C
C     **************************************************************
C
C     Modified to call BLAS.
C     Peter Mayes, NAG Central Office, October 1987.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPSMCH
      INTEGER           IFLAG, IGRADE, IV, LH, LIW, LJ, LW, M, N, NS
C     .. Array Arguments ..
      DOUBLE PRECISION  FJAC(LJ,N), FVEC(M), P(N), PHESD(NS), PHESL(LH),
     *                  RHS(NS), V(IV,N), W(LW), X(N)
      INTEGER           IW(LIW)
C     .. Subroutine Arguments ..
      EXTERNAL          LSF, LSFJSA
C     .. Local Scalars ..
      DOUBLE PRECISION  H, HINV, SUM, VTP
      INTEGER           I, J, JM1, JP1, JR, K, KKS, KS, KT, L, LFJAC,
     *                  LFVEC, LK, LL, LN, LX
C     .. External Subroutines ..
      EXTERNAL          F06FBF, DGEMV
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
C
C     ALLOCATE ALL THE LOCAL WORKSPACE
C
      LFVEC = N + 1
      LFJAC = LFVEC + M
      LX = LFJAC + M*N
      H = SQRT(EPSMCH)
      HINV = 1.0D+0/H
      KKS = 0
      KT = 1
C
C     FORM THE VECTOR TX = X + A SMALL FINITE DIFFERENCE
C
      DO 180 J = 1, NS
         JR = IGRADE + J
         LL = LX
         DO 20 I = 1, N
            W(LL) = X(I) + H*V(I,JR)
            LL = LL + 1
   20    CONTINUE
         IFLAG = 1
         CALL LSFJSA(IFLAG,M,N,LSF,W(LX),W(LFVEC),W(LFJAC),M,IW,LIW,W,
     *               LW)
         IF (IFLAG.LT.0) RETURN
         CALL F06FBF(NS,0.0D0,W,1)
         LL = LFJAC
         DO 60 K = 1, N
            DO 40 L = 1, M
               W(LL) = (W(LL)-FJAC(L,K))*HINV
               LL = LL + 1
   40       CONTINUE
   60    CONTINUE
         VTP = 0.0D+0
         LK = LFJAC
         LL = LK
         DO 100 L = 1, M
            CALL DGEMV('Transpose',N,NS,FVEC(L),V(1,IGRADE+1),IV,W(LL),
     *                 M,1.0D0,W,1)
            LL = LL + 1
            SUM = 0.0D+0
            LN = LK
            DO 80 K = 1, N
               SUM = SUM + P(K)*W(LN)
               LN = LN + M
   80       CONTINUE
            LK = LK + 1
C
C           THE MATRIX VTHP IS STORED IN RHS
C
            VTP = VTP + SUM*FVEC(L)
  100    CONTINUE
         RHS(J) = VTP
C
C        THE DIAGONAL ELEMENTS OF VTHV ARE STORED IN PHESD
C
         PHESD(J) = W(J)
         IF (J.EQ.NS) GO TO 140
         KS = KKS + J
         KKS = KS
         JP1 = J + 1
C
C        THE OFF-DIAGONAL ELEMENTS OF THE SYMMETRIC MATRIX VTHV
C        ARE STORED IN PHESL
C
         DO 120 I = JP1, NS
            PHESL(KS) = W(I)
            KS = KS + I - 1
  120    CONTINUE
  140    IF (J.EQ.1) GO TO 180
         JM1 = J - 1
C
C        SYMMETRIZE THE ELEMENTS OF VTHV ALREADY STORED IN PHESL
C
         DO 160 I = 1, JM1
            PHESL(KT) = (PHESL(KT)+W(I))*5.0D-1
            KT = KT + 1
  160    CONTINUE
  180 CONTINUE
      RETURN
C
C     END OF E04HEW   (FDHESS)
C
      END
