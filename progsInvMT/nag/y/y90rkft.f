      SUBROUTINE Y90RKF(TYPE,DTYPE,N,NBLOCK,BLOCK,A,IA,D,DIAG,COND,
     *                  SCALE,DETMAN,DETEXP,DIST,SEED,IWK1,IWK2,RWK1,
     *                  IRWK1,RWK2,IRWK2,RWK3,IRWK3,RWK4,IRWK4)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         =======================================
C         *  Y90RKF :  Random Matrix Generator  *
C         =======================================
C
C
C              ARGUMENTS
C              =========
C
C     TYPE   :  Integer, input.
C               Defines the type of diagonal blocks:
C               1  ==>  orthogonal matrices.
C               2  ==>  general matrix (unsymmetric).
C               3  ==>  symmetric matrix.
C
C     DTYPE  :  Defines how the vector D (used also to generate via
C               othogonal transformation random matrices) is generated.
C               0    ==> given as input
C               +/-1 ==> the elements of D have random values between
C                        DIAG(1) and DIAG(2).
C               +/-2 ==> use the parameter COND to generate elements
C                        using  D(I) = COND**(-(I-1)/(N-1)).
C               +/-3 ==> as above using (1-(I-1)/(N+1)+(I-1)/(COND*(N-1)
C               If DTYPE < 0 the order of the entries of D is reversed.
C
C     N      :  Integer, input.
C               Order of matrix to be generated.
C
C     NBLOCK :  Integer, output.
C               Number of blocks along the diagonal.
C
C     BLOCK  :  Integer array of DIMENSION (3,*), output.
C               For the i-th block BLOCK contains:
C               BLOCK(1,i)  ==>  Number of rows in the block.
C               BLOCK(2,i)  ==>  Number of columns in the block.
C               BLOCK(3,i)  ==>  Number of columns of overlap with the
C                                (i+1)-th block.
C
C     A      :  Real array of DIMENSION (IA,*), output.
C               Contains the matrix to be generated (square format).
C
C     IA     :  Integer, input.
C               Leading dimension of array A.  It must be  IA >= N.
C
C     D      :  Real array of DIMENSION (*), output.
C               Elements of the diagonal matrix (singular values) used
C               to generate the ABD matrix.
C
C     DIAG   :  Real array of DIMENSION (*), input.
C               DIAG(1) and DIAG(2) give the bounds on the random
C               elements of D if DTYPE = +/- 1 (see above).
C
C     COND   :  Real, input.
C               ABD matrix condition number.  Used if DTYPE is +/- 2 or
C               +/- 3 (see baove).
C
C     SCALE  :  Real, input.
C               Scaling factor.  All matrix elements (via its singular
C               values stored in D) are multiplied by SCALE.
C
C     DETMAN :  Real, output.
C               To avoid overflow the determinant of a square matrix
C               is calculated as: DETMAN*2.0**DETEXP, such that
C               0.0625 <= |DETMAN| < 1.0.
C
C     DETEXP :  Integer, output.
C               Exponent used to calculate the determinant (see DETMAN).
C
C     DIST   :  Integer, input.
C               Type of distribution used to generate random numbers.
C               1  ==>  Uniform (0,1)
C               2  ==>  Uniform (-1,1)
C               3  ==>  Normal (0,1)
C
C     SEED   :  Integer array of DIMENSION (4), Input/output.
C               Seeds for the random number generator.
C
C     IWK1, IWK2 :
C               Integer arrays of DIMENSION (*), used as workspace.
C
C     RWK1, RWK2, RWK3, RWK4 :
C               Real arrays of DIMENSION (IRWK1,*), (IRWK2,*),
C               (IRWK3,*), (IRWK4,*), respectively.  Used as workspace.
C
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  COND, DETMAN, SCALE
      INTEGER           DETEXP, DIST, DTYPE, IA, IRWK1, IRWK2, IRWK3,
     *                  IRWK4, N, NBLOCK, TYPE
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), D(*), DIAG(*), RWK1(IRWK1,*),
     *                  RWK2(IRWK2,*), RWK3(IRWK3,*), RWK4(IRWK4,*)
      INTEGER           BLOCK(3,*), IWK1(*), IWK2(*), SEED(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  X1, X2, X3
      INTEGER           DTYPE1, I, IC, IC1, IC2, IR1, IR2, J
C     .. External Functions ..
      DOUBLE PRECISION  Y90TBF
      EXTERNAL          Y90TBF
C     .. External Subroutines ..
      EXTERNAL          F06EPF, F06QHF, Y90RCF, Y90RGF, Y90RLF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN, NINT, DBLE, SQRT
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Generate the block structure
C
C-----------------------------------------------------------------------
      CALL Y90RLF(N,NBLOCK,BLOCK,SEED)
C-----------------------------------------------------------------------
C
C     Generate a random block-diagonal matrix
C
C-----------------------------------------------------------------------
C
C     1. Generate the matrix singular values
C
      CALL Y90RGF('Determinant',DTYPE,N,D,DIAG,COND,SCALE,DETMAN,DETEXP,
     *            DIST,SEED)
C
C     2. Random permutation of the singular values
C
      DO 20 I = 1, N
         J = MAX(1,MIN(NINT(DBLE(N)*Y90TBF(1,SEED)+HALF),N))
         IF (J.NE.I) THEN
            X1 = D(I)
            D(I) = D(J)
            D(J) = X1
         END IF
   20 CONTINUE
C
C     3. Initialize the block diagonal matrix
C
      CALL F06QHF('General',N,N,ZERO,ZERO,A,IA)
C
C     4. Generate the diagonal blocks
C
      DTYPE1 = 0
      IR2 = 0
      X1 = ONE
      DO 40 I = 1, NBLOCK
         IR1 = IR2 + 1
         IR2 = IR2 + BLOCK(1,I)
         CALL Y90RCF(TYPE,DTYPE1,BLOCK(1,I),BLOCK(1,I),A(IR1,IR1),IA,
     *               D(IR1),DIAG,COND,X1,DETMAN,DETEXP,DIST,SEED,
     *               'No permutation','No permutation',IWK1,IWK2,RWK1,
     *               IRWK1,RWK2,IRWK2,RWK3,IRWK3,RWK4,IRWK4)
   40 CONTINUE
C-----------------------------------------------------------------------
C
C     Generate randomly the Almost-Block-Diagonal matrix
C
C-----------------------------------------------------------------------
C
C     1. Use Givens' rotations of columns to obtain the matrix of
C        required structure starting from the bock-diagonal one.
C
      IR2 = 0
      DO 100 I = 1, NBLOCK
         IR1 = IR2 + 1
         IR2 = IR2 + BLOCK(1,I)
         IF (I.LE.1) THEN
            IC1 = 1
            IC2 = BLOCK(2,1)
         ELSE
            IC1 = IC2 - BLOCK(3,I-1) + 1
            IC2 = IC1 + BLOCK(2,I) - 1
            DO 60 IC = IR1 - 1, IC1, -1
               X1 = Y90TBF(2,SEED)
               X2 = Y90TBF(2,SEED)
               X3 = SQRT(X1*X1+X2*X2)
               X1 = X1/X3
               X2 = X2/X3
               CALL F06EPF(N,A(1,IC),1,A(1,IR1),1,X1,X2)
   60       CONTINUE
            DO 80 IC = IR1 + 1, IC1 + BLOCK(3,I-1) - 1
               X1 = Y90TBF(2,SEED)
               X2 = Y90TBF(2,SEED)
               X3 = SQRT(X1*X1+X2*X2)
               X1 = X1/X3
               X2 = X2/X3
               CALL F06EPF(N,A(1,IC),1,A(1,IR1),1,X1,X2)
   80       CONTINUE
         END IF
  100 CONTINUE
C-----------------------------------------------------------------------
C
C     End of Y90RKF
C
C-----------------------------------------------------------------------
      RETURN
      END
