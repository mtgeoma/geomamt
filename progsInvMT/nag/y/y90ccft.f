      SUBROUTINE Y90CCF(TYPE,DTYPE,M,N,A,IA,D,DIAG,COND,SCALE,DETMAN,
     *                  DETEXP,DIST,SEED,PVROW,PVCOL,IPVROW,IPVCOL,U,IU,
     *                  V,IV,WORK1,IWORK1,WORK2,IWORK2)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         =======================================
C         *  Y90CCF :  Random Matrix Generator  *
C         =======================================
C
C
C              ARGUMENTS
C              =========
C
C     TYPE   :  Type of random matrix to be generated.
C               0 ==>  only the vector D is generated
C               1 ==>  right or left orthogonal matrix generated
C               2 ==>  random rectangular matrix :  (Q**H)*D*Q
C               3 ==>  random hermitian positive definite square matrix
C               4 ==>  random symmetric positive definite square matrix
C
C     DTYPE  :  Defines how the vector D (used also to generate via
C               orthogonal transformation random matrices) is generated.
C               0    ==> given as input
C               +/-1 ==> the elements of D have values with real parts
C                        between the real parts of DIAG(1),DIAG(2),
C                        likewise for the imaginary parts
C               +/-2 ==> use the parameter COND to generate real
C                        elements using  D(I) = COND**(-(I-1)/(N-1))
C               +/-3 ==> as above using (1-(I-1)/(N+1)+(I-1)/(COND*(N-1)
C               +/-4 ==> as for +/-1 but real values are multiplied by
C                        random complex numbers on the unit circle
C               +/-5 ==> as for +/-2 but real values are multiplied by
C                        random complex numbers on the unit circle
C               If DTYPE < 0 the order of the entries of D is reversed.
C
C     M      :  Number of rows in the matrix to be generated.
C
C     N      :  Number of columns in the matrix to be generated.
C
C     A      :  Array storing the generated random matrix if TYPE =/= 0.
C               If |TYPE| = 1 it stores a left orthogonal matrix for
C               TYPE > 0, a right orthogonal for TYPE < 0.
C
C     IA     :  First dimension of the array A .
C
C     D      :  If DTYPE=0 it stores the diagonal matrix D on entry.
C               Otherwise it stores the generated diagonal elements of D
C               for random matrix generation (TYPE=/=0) or a random
C               vector (|DTYPE|=1 and TYPE=0).
C
C     DIAG   :  Contains on entry limiting values for the elements of D.
C
C     COND   :  Condition number: used to generate the singular values
C               stored in D (|DTYPE|=2,3,4 or 5).
C
C     SCALE  :  Scaling factor which multiplies all elements of D.
C
C     DETMAN :  To avoid overflow the determinant of a square matrix
C               is calculated as: DETMAN*2.0**DETEXP, such that
C               0.0625 <= |DETMAN| < 1.0.
C
C     DETEXP :  Exponent used to calculate the determinant (see DETMAN).
C
C     DIST   :  Type of distribution used to generate random numbers.
C               For real random numbers (|DTYPE| = 1):
C               1 ==>  Uniform (0,1)
C               2 ==>  Uniform (-1,1)
C               3 ==>  Normal (0,1)
C               For complex random numbers (|DTYPE| = 4 or 5):
C               1 ==>  real and imaginary parts uniform (0,1)
C               2 ==>      "        "          "        (-1,1)
C               3 ==>  random comple in the disc |z| < 1
C               4 ==>  real and imaginary parts normal (0,1)
C               5 ==>  random complex on the circle |z| = 1
C
C     SEED   :  Seeds for the random number generator.
C
C     PVROW  :  Switch for permutation of rows:
C               'P' ==>  permute using IPVROW.
C
C     PVCOL  :  Switch for permutation of columns:
C               'P' ==>  permute using IPVCOL.
C
C     IPVROW :  Specifies the permutation of rows used.  It is the
C               inverseof LINPACK pivoting.
C
C     IPVCOL :  Specifies the permutation of columns used.  It is the
C               inverseof LINPACK pivoting.
C
C     U      :  Working space.
C
C     V      :  Working space.
C
C     WORK1  :  Working space
C
C     WORK2  :  Working space.
C
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      COMPLEX*16        CZERO, CONE
      PARAMETER         (CZERO=(0.0D0,0.0D0),CONE=(1.0D0,0.0D0))
C     .. Scalar Arguments ..
      COMPLEX*16        DETMAN, SCALE
      DOUBLE PRECISION  COND
      INTEGER           DETEXP, DIST, DTYPE, IA, IU, IV, IWORK1, IWORK2,
     *                  M, N, TYPE
      CHARACTER*1       PVCOL, PVROW
C     .. Array Arguments ..
      COMPLEX*16        A(IA,*), D(*), DIAG(*), U(IU,*), V(IV,*),
     *                  WORK1(IWORK1,*), WORK2(IWORK2,*)
      INTEGER           IPVCOL(*), IPVROW(*), SEED(4)
C     .. Local Scalars ..
      INTEGER           I, J, ND
      CHARACTER*1       DETREQ
C     .. External Functions ..
      LOGICAL           Y90WAF
      EXTERNAL          Y90WAF
C     .. External Subroutines ..
      EXTERNAL          F06GDF, F06GGF, F06ZAF, Y90CDF, Y90CGF, Y90DHF,
     *                  Y90DKF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Initialize
C
C-----------------------------------------------------------------------
      ND = MIN(M,N)
C-----------------------------------------------------------------------
C
C     Generate diagonal
C
C-----------------------------------------------------------------------
      IF (TYPE.NE.1 .AND. DTYPE.NE.0) THEN
         IF (M.EQ.N) THEN
            DETREQ = 'D'
         ELSE
            DETREQ = 'N'
         END IF
         CALL Y90CGF(DETREQ,DTYPE,N,D,DIAG,COND,SCALE,DETMAN,DETEXP,
     *               DIST,SEED)
      END IF
C-----------------------------------------------------------------------
C
C     Generate Random or Orthogonal Matrix
C
C-----------------------------------------------------------------------
      IF (TYPE.NE.0) THEN
C
         IF (M.GT.1) THEN
            DO 40 I = 1, M
               DO 20 J = 1, N
                  A(I,J) = CZERO
   20          CONTINUE
   40       CONTINUE
         END IF
C
C     Orthogonal matrix
C
         IF (ABS(TYPE).EQ.1) THEN
            IF (TYPE.GE.0) THEN
               CALL Y90CDF('Left','Initialize',M,N,A,IA,WORK1,WORK2,
     *                     SEED)
            ELSE
               CALL Y90CDF('Right','Initialize',M,N,A,IA,WORK1,WORK2,
     *                     SEED)
            END IF
C
C     General random matix
C
         ELSE IF (TYPE.EQ.2) THEN
            DO 60 I = 1, ND
               A(I,I) = D(I)
   60       CONTINUE
            CALL Y90CDF('Left','No Initialize',M,N,A,IA,WORK1,WORK2,
     *                  SEED)
            CALL Y90CDF('Right','No Initialize',M,N,A,IA,WORK1,WORK2,
     *                  SEED)
C
C     Hermitian matrix
C
         ELSE IF (TYPE.EQ.3) THEN
            CALL Y90CDF('Left','Initialize',M,N,U,IU,WORK1,WORK2,SEED)
            CALL Y90DHF('C',M,N,U,IU,V,IV)
            DO 80 I = 1, M
               CALL F06GDF(M,D(I),U(1,I),1)
   80       CONTINUE
            CALL F06ZAF('N','N',M,M,M,CONE,U,IU,V,IV,CZERO,A,IA)
C
C     Symmetrize the matrix A
C
            CALL Y90DKF('U',A,IA,N)
C
C     Symmetric Complex matrix
C
         ELSE IF (TYPE.EQ.4) THEN
C
            CALL Y90CDF('Left','Initialize',M,N,U,IU,WORK1,WORK2,SEED)
            CALL Y90DHF('T',M,N,U,IU,V,IV)
            DO 100 I = 1, M
               CALL F06GDF(M,D(I),U(1,I),1)
  100       CONTINUE
            CALL F06ZAF('N','N',M,M,M,CONE,U,IU,V,IV,CZERO,A,IA)
C
C     Symmetrize the matrix A (do not use Y90DKF because it
C     sets imaginary parts of diagonal elements to zero)
C
            DO 140 J = 1, N
               DO 120 I = J + 1, N
                  A(I,J) = A(J,I)
  120          CONTINUE
  140       CONTINUE
C
         END IF
C-----------------------------------------------------------------------
C
C     Permute Rows and/or columns when required
C
C-----------------------------------------------------------------------
         IF (Y90WAF(PVROW,'P')) THEN
            DO 160 I = M, 1, -1
               CALL F06GGF(N,A(I,1),IA,A(IPVROW(I),1),IA)
  160       CONTINUE
         END IF
C
         IF (Y90WAF(PVCOL,'P')) THEN
            DO 180 I = N, 1, -1
               CALL F06GGF(N,A(1,I),1,A(1,IPVCOL(I)),1)
  180       CONTINUE
         END IF
C
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90CCF
C
C-----------------------------------------------------------------------
      RETURN
      END
