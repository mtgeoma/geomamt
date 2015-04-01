      SUBROUTINE Y90RRF(N,NZ,A,ICN,IRN,OCCUP,NORMA,D,HIT,LDHIT,IW,AMAT,
     *                  LDAMAT,SEED)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C-----------------------------------------------------------------------
C
C         =======================================
C         *  Y90RRF :  Sparse Matrix Generator  *
C         =======================================
C
C
C              DESCRIPTION
C              ===========
C
C     Y90RRF generates a random real symmetric sparse matrix of known
C     eigenvalues (singular values) using a sequence of random Givens
C     rotation to obtain the required level of off-diagonal occupancy.
C
C
C              ARGUMENTS
C              =========
C
C     N      :  Integer, input.
C               Order of the matrix.
C     NZ     :  Integer, output.
C               Number of non-zeros in the generated matrix.
C     A      :  Real, output,
C               Arrays of dimension (*).
C               Non-zero elements of the generated matrix A (by column).
C     ICN    :  Integer, output.
C               Arrays of dimension (*).
C               Column indices of the elements stored in A.
C     IRN    :  Integer, output.
C               Arrays of dimension (*).
C               Row indices of the elements stored in A.
C     OCCUP  :  Real, input.
C               Occupancy level (in percent) of the off-diagonal
C               elements required, where 0 is no off-diagonal non-zero
C               elements, 100 all off-diagonal non-zero elements.
C     NORMA  :  Real, output.
C               1-norm of the generated matrix.
C     D      :  Real, input.
C               Array of dimension (*).
C               D contains the singular (eigen) values of the matrix A.
C     HIT    :  Logical, work space.
C               Two-dimensional array of dimension (LDHIT,*)
C     LDHIT  :  Integer, input.
C               Leading dimension of the array HIT.  LDHIT >= N.
C     IW     :  Integer, work space.
C               Array of dimension (N,*).
C     AMAT   :  Real, output,
C               Two-dimensional arrays of dimension (LDAMAT,*).
C               The generated sparse matrix stored as a dense matrix.
C     LDAMAT :  Integer, input.
C               Leading dimension of the array AMAT.  LDAMAT >= N.
C     SEED   :  Integer, input/output.
C               Seeds of the random number generator.  Updated by each
C               call to it.
C
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF
      PARAMETER         (ZERO=0.0D0,HALF=0.5D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  NORMA, OCCUP
      INTEGER           LDAMAT, LDHIT, N, NZ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*), AMAT(LDAMAT,*), D(*)
      INTEGER           ICN(*), IRN(*), IW(N,*), SEED(4)
      LOGICAL           HIT(LDHIT,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  CTHETA, PI, STHETA, TEMP, THETA
      INTEGER           I, IC1, IC2, IN, IROT, IZ, J, K, NK, NZMIN, OCC
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF, Y90TBF
      EXTERNAL          X01AAF, Y90TBF
C     .. External Subroutines ..
      EXTERNAL          F06DBF, F06EFF, F06EPF, F06QHF, Y90RRX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, COS, DBLE, MAX, MIN, NINT, SIN
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Initialize the matrix A and other items
C
C-----------------------------------------------------------------------
      CALL F06QHF('G',N,N,ZERO,ZERO,AMAT,LDAMAT)
      CALL F06EFF(N,D,1,AMAT,LDAMAT+1)
C
      OCC = 0
C
      DO 40 J = 1, N
         DO 20 I = 1, N
            HIT(I,J) = .FALSE.
   20    CONTINUE
         HIT(J,J) = .TRUE.
   40 CONTINUE
C
      PI = X01AAF(PI)
C-----------------------------------------------------------------------
C
C     Loop over random rotations
C
C-----------------------------------------------------------------------
      DO 180 IROT = 1, 10000
C-----------------------------------------------------------------------
C
C        Count the number of non-zeros in each column
C
C-----------------------------------------------------------------------
         CALL F06DBF(N,0,IW(1,1),1)
         CALL F06DBF(N,0,IW(1,2),1)
         CALL F06DBF(N,0,IW(1,3),1)
C
         NZMIN = N + 1
         DO 80 J = 1, N
            NZ = 0
            DO 60 I = 1, N
               IF (AMAT(I,J).NE.ZERO) NZ = NZ + 1
   60       CONTINUE
            NZMIN = MIN(NZMIN,NZ)
            IW(NZ,3) = IW(NZ,3) + 1
            IN = IW(NZ,1)
            IW(NZ,1) = J
            IW(J,2) = IN
   80    CONTINUE
C-----------------------------------------------------------------------
C
C        Select randomly the first column from the columns of
C        minimum degree
C
C-----------------------------------------------------------------------
         NK = IW(NZMIN,3)
         I = MIN(NK,MAX(1,NINT(DBLE(NK)*Y90TBF(1,SEED)+HALF)))
         IC1 = IW(NZMIN,1)
         DO 100 J = 1, I - 1
            IC1 = IW(IC1,2)
  100    CONTINUE
C-----------------------------------------------------------------------
C
C        Now select the other column
C
C-----------------------------------------------------------------------
         IC2 = 0
         DO 120 IZ = NZMIN, N
            IF (IW(IZ,3).GT.0) THEN
               CALL Y90RRX(IC1,IZ,IW,IW(1,2),IW(1,3),IW(1,4),NK,HIT,
     *                     LDHIT)
               IF (NK.GT.0) THEN
                  I = MIN(NK,MAX(1,NINT(DBLE(NK)*Y90TBF(1,SEED)+HALF)))
                  IC2 = IW(I,4)
                  GO TO 140
               END IF
            END IF
  120    CONTINUE
  140    CONTINUE
C-----------------------------------------------------------------------
C
C     Count the number of new non-zero elements
C
C-----------------------------------------------------------------------
         IF (IC2.GT.0) THEN
C
            K = 0
            DO 160 I = 1, N
               IF (((AMAT(IC1,I).EQ.ZERO) .AND. (AMAT(IC2,I).NE.ZERO))
     *              .OR. ((AMAT(IC2,I).EQ.ZERO) .AND. (AMAT(IC1,I)
     *             .NE.ZERO))) K = K + 1
  160       CONTINUE
            IF (AMAT(IC1,IC2).EQ.ZERO) THEN
               K = K + K - 2
            ELSE
               K = K + K
            END IF
C-----------------------------------------------------------------------
C
C        Now apply the rotations
C
C-----------------------------------------------------------------------
            IF (100*DBLE(OCC+K)/DBLE(N*(N-1)).LE.OCCUP) THEN
C
               OCC = OCC + K
C
               HIT(IC1,IC2) = .TRUE.
               HIT(IC2,IC1) = .TRUE.
               THETA = PI*Y90TBF(2,SEED)
               CTHETA = COS(THETA)
               STHETA = SIN(THETA)
               CALL F06EPF(N,AMAT(IC1,1),LDAMAT,AMAT(IC2,1),LDAMAT,
     *                     CTHETA,STHETA)
               CALL F06EPF(N,AMAT(1,IC1),1,AMAT(1,IC2),1,CTHETA,STHETA)
C
            ELSE
C
               GO TO 200
C
            END IF
C
         END IF
C-----------------------------------------------------------------------
C
C     End of loop over random rotations
C
C-----------------------------------------------------------------------
  180 CONTINUE
  200 CONTINUE
C-----------------------------------------------------------------------
C
C     Collect the non-zeros in the upper triangle and compute the norm
C
C-----------------------------------------------------------------------
      NZ = 0
      NORMA = ZERO
      DO 260 J = 1, N
         TEMP = ZERO
         DO 220 I = 1, J
            IF (AMAT(I,J).NE.ZERO) THEN
               NZ = NZ + 1
               A(NZ) = AMAT(I,J)
               ICN(NZ) = J
               IRN(NZ) = I
               TEMP = TEMP + ABS(A(NZ))
            END IF
  220    CONTINUE
         DO 240 I = J + 1, N
            IF (AMAT(I,J).NE.ZERO) TEMP = TEMP + ABS(AMAT(I,J))
  240    CONTINUE
         NORMA = MAX(NORMA,TEMP)
  260 CONTINUE
C-----------------------------------------------------------------------
C
C     End of Y90RRF
C
C-----------------------------------------------------------------------
      RETURN
      END
