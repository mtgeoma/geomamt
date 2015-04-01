      SUBROUTINE Y90RMF(DET,PERM,BSTRUC,BDIAG,TYPE,M,N,DIST,DTYPE,COND,
     *                  SCALE,DETMAN,DETEXP,D,DIAG,A,NA,ICOL,IROW,V,W,
     *                  NBLOCK,BLOCK,SEED)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ============================================
C         *  Y90RMF Generate random sparse matrices  *
C         ============================================
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, HALF
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  COND, DETMAN, SCALE
      INTEGER           DETEXP, DIST, DTYPE, M, N, NA, NBLOCK, TYPE
      CHARACTER*1       BDIAG, BSTRUC, DET, PERM
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*), D(*), DIAG(*), V(*), W(*)
      INTEGER           BLOCK(*), ICOL(*), IROW(*), SEED(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  COST, PI, SINT, TEMP, THETA
      INTEGER           BLOCK1, BLOCK2, I, IB, IB1, IB2, II, IJ, J, K,
     *                  K1, K2, KT, L, L1, L2, MNMAX, MNMIN, NB, NI, NK,
     *                  NNVSTR, NROT1, NROT2, NTOT, NV
C     .. External Functions ..
      DOUBLE PRECISION  F06EJF, X01AAF, Y90TBF
      LOGICAL           Y90WAF
      EXTERNAL          F06EJF, X01AAF, Y90TBF, Y90WAF
C     .. External Subroutines ..
      EXTERNAL          F06EFF, F06EGF, F06EPF, F06FBF, F06PAF, F06PMF,
     *                  Y90RGF, Y90RNF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, COS, MAX, MIN, NINT, DBLE, SIGN, SIN
C     .. Statement Functions ..
      INTEGER           INTRND
C     .. Statement Function definitions ..
      INTRND(L1,L2) = MIN(L2,MAX(L1,NINT(DBLE(L2-L1+1)*Y90TBF(1,SEED)
     *                +HALF)+L1-1))
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     1. Generate diagonal matrix
C
C-----------------------------------------------------------------------
      MNMAX = MAX(M,N)
      MNMIN = MIN(M,N)
C
      IF (DTYPE.NE.0) THEN
         CALL Y90RGF(DET,DTYPE,MNMIN,D,DIAG,COND,SCALE,DETMAN,DETEXP,
     *               DIST,SEED)
         DO 20 J = 1, MNMIN
            IROW(J) = 0
            K = MAX(1,MIN(MNMIN,NINT(DBLE(MNMIN)*Y90TBF(1,SEED)+HALF)))
            IF (K.NE.J) THEN
               TEMP = D(J)
               D(J) = D(K)
               D(K) = TEMP
            END IF
   20    CONTINUE
      END IF
C-----------------------------------------------------------------------
C
C     2. Apply random Givens rotations to form a first sparse matrix
C
C        If a symmetric positive-definite sparse matrix was required,
C        the process terminates here.
C
C        Otherwise subsequent steps are carried out
C
C-----------------------------------------------------------------------
      BLOCK1 = MIN(3,(N+5)/6)
      BLOCK2 = MAX(4,(N+2)/3)
      IF (TYPE.LE.1) THEN
         IF (Y90WAF(BDIAG,'Diagonal')) THEN
            NROT1 = 0
            NROT2 = 2*MNMAX/3
         ELSE
            NROT1 = MNMAX/2
            NROT2 = MNMAX/2
         END IF
      ELSE
         IF (Y90WAF(BDIAG,'Diagonal')) THEN
            NROT1 = 0
            NROT2 = 2*MNMAX/3
         ELSE
            NROT1 = MNMAX
            NROT2 = 0
         END IF
      END IF
C
      CALL F06FBF(M*N,ZERO,A,1)
      CALL F06EFF(MNMIN,D,1,A,M+1)
C
      PI = X01AAF(PI)
C
      IF ( .NOT. Y90WAF(BDIAG,'Diagonal')) THEN
         IF (TYPE.LE.1) THEN
            NK = MNMIN - 1
            DO 40 J = MNMIN + 1, MNMAX
               THETA = PI*Y90TBF(2,SEED)
               COST = COS(THETA)
               SINT = SIN(THETA)
               NK = NK + 1
               K = MAX(1,MIN(NK,NINT(DBLE(NK)*Y90TBF(1,SEED)+HALF)))
               IF (M.GE.N) THEN
                  CALL F06EPF(N,A(K),M,A(J),M,COST,SINT)
               ELSE
                  CALL F06EPF(M,A(M*(K-1)+1),1,A(M*(J-1)+1),1,COST,SINT)
               END IF
   40       CONTINUE
         END IF
C
         NV = 0
         DO 100 J = 1, NROT1
            THETA = PI*Y90TBF(2,SEED)
            COST = COS(THETA)
            SINT = SIN(THETA)
C
   60       CONTINUE
            K = INTRND(1,M)
            IF ((IROW(K).GE.1) .AND. NV.LE.M/2) GO TO 60
C
            IF (IROW(K).LE.0) NV = NV + 1
            IROW(K) = 1
C
   80       CONTINUE
            L = INTRND(1,M)
            IF ((L.EQ.K) .OR. ((IROW(L).GE.1) .AND. NV.LE.M/2)) GO TO 80
C
            NV = NV + 1
            IROW(L) = 1
            IF (TYPE.LE.1) THEN
               CALL F06EPF(N,A(K),M,A(L),M,COST,SINT)
            ELSE
               CALL F06EPF(N,A(K),M,A(L),M,COST,SINT)
               CALL F06EPF(M,A(M*(K-1)+1),1,A(M*(L-1)+1),1,COST,SINT)
            END IF
  100    CONTINUE
C-----------------------------------------------------------------------
C
C     Form now a block triangular sparse matrix (unsymmetric sparse
C     matrices only)
C
C-----------------------------------------------------------------------
C
C     1. Apply Householder's reflections to triangularise matrix
C
         IF (TYPE.LE.1) THEN
            DO 120 I = 1, MNMIN
               NI = N - I + 1
               II = M*(I-1) + I
               CALL F06EFF(NI,A(II),M,V,1)
               TEMP = F06EJF(NI,V,1)
               V(1) = V(1) + SIGN(TEMP,V(1))
               TEMP = -ONE/(TEMP*(TEMP+ABS(A(II))))
C
               CALL F06PAF('N',M,NI,ONE,A(M*(I-1)+1),M,V,1,ZERO,W,1)
               CALL F06PMF(M,NI,TEMP,W,1,V,1,A(M*(I-1)+1),M)
               CALL F06FBF(NI-1,ZERO,A(II+M),M)
  120       CONTINUE
         END IF
      END IF
C
C     2. Define block structure
C
C        a. Generate block sizes
C
      IF (Y90WAF(BSTRUC,'Block')) THEN
         NTOT = 0
         DO 140 I = 1, MNMIN
            NBLOCK = NINT(DBLE(BLOCK2-BLOCK1+1)*Y90TBF(1,SEED)-0.5D0) +
     *               BLOCK1
            NBLOCK = MIN(BLOCK2,MAX(BLOCK1,NBLOCK))
            IF (NTOT+NBLOCK.LT.MNMIN) THEN
               NTOT = NTOT + NBLOCK
               BLOCK(I) = NBLOCK
            ELSE
               BLOCK(I) = MNMIN - NTOT
               NBLOCK = I
               GO TO 160
            END IF
  140    CONTINUE
  160    CONTINUE
         DO 180 J = 1, NBLOCK
            K = INTRND(1,NBLOCK)
            IF (K.NE.J) THEN
               L = BLOCK(J)
               BLOCK(J) = BLOCK(K)
               BLOCK(K) = L
            END IF
  180    CONTINUE
      END IF
C
C        b. Apply within-block Givens' rotations
C
      DO 280 J = 1, NROT2
  200    CONTINUE
         K = INTRND(1,MNMIN)
         IB2 = 0
         DO 220 IB = 1, NBLOCK
            IB1 = IB2 + 1
            IB2 = IB2 + BLOCK(IB)
            IF ((K.GE.IB1) .AND. (K.LE.IB2)) GO TO 240
  220    CONTINUE
  240    CONTINUE
         IF (IB2-IB1.LE.0) GO TO 200
C
  260    CONTINUE
         L = INTRND(IB1,IB2)
         IF (K.EQ.L) GO TO 260
         THETA = PI*Y90TBF(2,SEED)
         COST = COS(THETA)
         SINT = SIN(THETA)
         IF (TYPE.LE.1) THEN
            IF (Y90TBF(2,SEED).LE.ZERO) THEN
               CALL F06EPF(IB2,A(K),M,A(L),M,COST,SINT)
            ELSE
               CALL F06EPF(M-IB1+1,A(M*(K-1)+IB1),1,A(M*(L-1)+IB1),1,
     *                     COST,SINT)
            END IF
         ELSE
            CALL F06EPF(IB2,A(K),N,A(L),N,COST,SINT)
            CALL F06EPF(N-IB1+1,A(N*(K-1)+IB1),1,A(N*(L-1)+IB1),1,COST,
     *                  SINT)
         END IF
C
  280 CONTINUE
C-----------------------------------------------------------------------
C
C        c. Make sure all blocks are irreducible
C
C-----------------------------------------------------------------------
      K2 = 0
      DO 320 I = 1, NBLOCK
         K1 = K2 + 1
         NB = BLOCK(I)
         K2 = K2 + NB
C
         CALL Y90RNF(NB,A(M*(K1-1)+K1),M,NNVSTR,IROW,IROW(NB+2),ICOL,
     *               ICOL(NB+1),ICOL(2*NB+1),ICOL(3*NB+1),ICOL(4*NB+1),
     *               ICOL(5*NB+1),ICOL(6*NB+1))
C
         KT = NNVSTR
         DO 300 J = NNVSTR - 1, 1, -1
            K = IROW(NB+1+INTRND(IROW(J),IROW(J+1)-1)) + K1 - 1
            L = IROW(NB+1+INTRND(IROW(KT),IROW(KT+1)-1)) + K1 - 1
            IF (TYPE.LE.1) THEN
               IF (Y90TBF(1,SEED).LE.0.0D0) THEN
                  THETA = PI*Y90TBF(2,SEED)
                  COST = COS(THETA)
                  SINT = SIN(THETA)
                  CALL F06EPF(K2,A(K),M,A(L),M,COST,SINT)
               ELSE
                  THETA = PI*Y90TBF(2,SEED)
                  COST = COS(THETA)
                  SINT = SIN(THETA)
                  CALL F06EPF(M-K1+1,A(M*(K-1)+K1),1,A(M*(L-1)+K1),1,
     *                        COST,SINT)
               END IF
            ELSE
               THETA = PI*Y90TBF(2,SEED)
               COST = COS(THETA)
               SINT = SIN(THETA)
               CALL F06EPF(N,A(K),M,A(L),M,COST,SINT)
               CALL F06EPF(M,A(M*(K-1)+1),1,A(M*(L-1)+1),1,COST,SINT)
            END IF
            KT = J
  300    CONTINUE
  320 CONTINUE
C-----------------------------------------------------------------------
C
C     4. Random permutation of rows and columns
C
C-----------------------------------------------------------------------
      IF (Y90WAF(PERM,'Permutation')) THEN
         IF (TYPE.LE.1) THEN
            DO 340 J = 1, N
               K = MAX(1,MIN(N,NINT(DBLE(N)*Y90TBF(1,SEED)+HALF)))
               IF (K.NE.J) CALL F06EGF(M,A(M*(J-1)+1),1,A(M*(K-1)+1),1)
  340       CONTINUE
            DO 360 J = 1, M
               K = MAX(1,MIN(M,NINT(DBLE(M)*Y90TBF(1,SEED)+HALF)))
               IF (K.NE.J) CALL F06EGF(N,A(J),M,A(K),M)
  360       CONTINUE
         ELSE
            DO 380 J = 1, N
               K = MAX(1,MIN(N,NINT(DBLE(N)*Y90TBF(1,SEED)+HALF)))
               IF (K.NE.J) THEN
                  CALL F06EGF(M,A(M*(J-1)+1),1,A(M*(K-1)+1),1)
                  CALL F06EGF(N,A(J),M,A(K),M)
               END IF
  380       CONTINUE
         END IF
      END IF
C-----------------------------------------------------------------------
C
C     5. Copy sparse matrix into final format
C
C-----------------------------------------------------------------------
      NA = 0
      IJ = 0
      DO 420 J = 1, N
         DO 400 I = 1, M
            IJ = IJ + 1
            IF (A(IJ).NE.ZERO) THEN
               NA = NA + 1
               A(NA) = A(IJ)
               ICOL(NA) = J
               IROW(NA) = I
            END IF
  400    CONTINUE
  420 CONTINUE
C-----------------------------------------------------------------------
C
C     End of Y90RMF
C
C-----------------------------------------------------------------------
      RETURN
      END
