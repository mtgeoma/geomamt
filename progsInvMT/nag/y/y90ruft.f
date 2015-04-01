      SUBROUTINE Y90RUF(DTYPE,N,DENS,NNZ,A,IROW,ICOL,IDIMA,LAMBDA,
     *                  LAMBND,COND,SCALE,DIST,SEED,IWORK,WORK)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
*-----------------------------------------------------------------------
*
*         ========================================================
*         *  Y90RUF :  Random Sparse Symmetric Matrix Generator  *
*         ========================================================
*
*
*     ARGUMENT LIST
*     =============
*
*  DTYPE     -  Integer (input)
*               Flag used for the generation of the eigenvalues of the
*               matrix A:
*               0 ==> these must be supplied in input in D
*               1 ==> SCALE*random(LAMBND(1),LAMBND(2)) (see DIST below)
*               2 ==> LAMBDA(k) = SCALE*COND**(-(k-1)/(n-1))
*               3 ==> LAMBDA(k) = SCALE*
*                                 ((1-(k-1)/(n-1))+(k-1)/(k-1))/COND
*  N         -  Integer (input)
*               The order of the matrix
*  DENS      -  Real (input)
*               Density to be achieved in percentage
*               (0.0 <= DENS <= 100.0)
*  NNZ       -  Integer (output)
*               The number of non-zeros
*  A(*)      -  Real array of size = IDIMA (output)
*               The non-zero elements of the matrix A
*  IROW(*)   -  Integer array of size = IDIMA (output)
*               The row indices of the non-zero elements of the matrix A
*  ICOL(*)   -  Integer array of size = IDIMA(output)
*               The column indices of the non-zero elements of the
*               matrix A
*  IDIMA     -  Integer (input)
*               The size of the arrays A, IROW and ICOL.  IDIMA >= NNZ
*  LAMBDA(*) -  Real array of size >= N (input/output)
*               On entry: if DTYPE = 0, the eigenvalues
*               On exit:  if DTYPE = 0, unchanged,
*                         otherwise, the generated eigenvalues: in this
*                         case LAMBDA(k) > 0 for all k
*  LAMBND(*) -  Real array of size >= 2 (input)
*               LAMBND(1), LAMBND(2) contain the two values used in
*               the generation of D when DTYPE = 1
*               Unreferenced if DTYPE =/= 1
*  COND      -  Real (input)
*               The condition number required
*               Unreferenced if DTYPE =/= 2 or 3
*  SCALE     -  Real (input)
*               Scaling factor for the eigenvalues
*               Unreferenced if DTYPE = 0
*  DIST      -  Integer (input)
*               Type of random distribution for the eigenvalues
*               1, 2 => uniform (LAMBND(1),LAMBND(2))
*               3    => LAMBND(1) + LAMBND(2)*normal(1,0) ???
*               Unreferenced if DTYPE =/= 1
*  SEED(4)   -  Integer array of size = 4 (input/output)
*               Seeds for the random number generator
*  IWORK(*)  -  Integer workspace of size >= 3*N
*  WORK(*)   -  Real workspace of size >= 2*N
*
*-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ONE, HUND
      PARAMETER         (ONE=1.0D0,HUND=1.0D2)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  COND, DENS, SCALE
      INTEGER           DIST, DTYPE, IDIMA, N, NNZ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*), LAMBDA(*), LAMBND(2), WORK(*)
      INTEGER           ICOL(*), IROW(*), IWORK(*), SEED(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  DETMAN, X
      INTEGER           DETEXP, I, IDUP, IFAIL, J, K, L, MAXDUP, NDUP,
     *                  NK
      LOGICAL           DIAG
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF, Y90TBF
      EXTERNAL          X02AMF, Y90TBF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, F06DBF, F06EFF, F06FBF, M01CAF, Y90RAF,
     *                  Y90RUZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, INT, MAX, MIN
C     .. Executable Statements ..
*-----------------------------------------------------------------------
*
*     Generate a tridiagonal matrix
*
*-----------------------------------------------------------------------
      DIAG = .FALSE.
      NDUP = 0
      MAXDUP = 1
      IWORK(1) = 1
*
*     DTYPE = 0.  If there are repeated eigenvalues, generate a
*                 piecewise irreducible tridiagonal matrix
*
      IF (DTYPE.LE.0) THEN
         CALL F06EFF(N,LAMBDA,1,WORK,1)
         IFAIL = 0
         CALL M01CAF(WORK,1,N,'Ascending',IFAIL)
*
         NDUP = 0
         CALL F06DBF(N,0,IWORK,1)
         X = MIN(-WORK(1),-ONE)
         DO 20 I = 1, N
            IF (WORK(I).EQ.X) THEN
               IWORK(NDUP) = IWORK(NDUP) + 1
            ELSE
               NDUP = NDUP + 1
               IWORK(NDUP) = 1
               X = WORK(I)
            END IF
   20    CONTINUE
         IWORK(NDUP+1) = 0
         MAXDUP = IWORK(NDUP)
         DO 40 I = 1, NDUP - 1
            MAXDUP = MAX(MAXDUP,IWORK(I))
            IWORK(NDUP+I+1) = IWORK(NDUP+I) + IWORK(I)
   40    CONTINUE
         IF (MAXDUP.LE.1) CALL DCOPY(N,LAMBDA,1,WORK,1)
*
*        All identical eigenvalues: the matrix is diagonal
*
         IF (MAXDUP.EQ.N) THEN
            DIAG = .TRUE.
            X = LAMBDA(1)
*
*        The matrix is not diagonal
*
         ELSE
            IDUP = 1
            I = 0
            IWORK(2*NDUP+1) = 1
            DO 80 K = 1, MAXDUP
               DO 60 L = IDUP, NDUP
                  IF (IWORK(L).GT.0) THEN
                     I = I + 1
                     IWORK(L) = IWORK(L) - 1
                     IWORK(NDUP+L) = IWORK(NDUP+L) + 1
                     WORK(N+I) = WORK(IWORK(NDUP+L))
                     IF ((IWORK(L).LE.0) .AND. (IDUP+1.EQ.L)) IDUP = L
                  ELSE IF (IDUP+1.EQ.L) THEN
                     IDUP = L
                  END IF
   60          CONTINUE
               IF (K.LE.MAXDUP) IWORK(2*NDUP+K+1) = I + 1
   80       CONTINUE
         END IF
*
*     DTYPE = 1.  If LAMBND(1) = LAMBND(2), the matrix is diagonal
*                 If LAMBND(1) =/= LAMBND(2), the matrix is tridiagonal
*
      ELSE IF (DTYPE.EQ.1) THEN
         IF (LAMBND(1).EQ.LAMBND(2)) THEN
            DIAG = .TRUE.
            X = SCALE*LAMBND(1)
            CALL F06FBF(N,X,LAMBDA,1)
         END IF
*
*     DTYPE = 2.  If COND = 1.0, the matrix is diagonal
*                 If COND =/= 1.0, the matrix is tridiagonal
*
      ELSE IF (DTYPE.EQ.2) THEN
         IF (COND.EQ.ONE) THEN
            DIAG = .TRUE.
            X = SCALE
            CALL F06FBF(N,X,LAMBDA,1)
         END IF
      END IF
*
*     Generate a diagonal matrix
*
      IF (DIAG) THEN
         COND = ONE
         NNZ = N
         CALL F06FBF(N,X,A,1)
         DO 100 I = 1, N
            IROW(I) = I
            ICOL(I) = I
  100    CONTINUE
*
*     Generate a piece-wise tridiagonal matrix
*
      ELSE
         NNZ = 1
         I = 0
         DO 140 K = 1, MAXDUP
            L = IWORK(2*NDUP+K)
            IF (K.LT.MAXDUP) THEN
               NK = IWORK(2*NDUP+K+1) - L
            ELSE
               NK = N - L + 1
            END IF
            IF (NK.LE.1) THEN
               A(NNZ) = WORK(N+IWORK(2*NDUP+K))
            ELSE
               CALL Y90RAF('Symmetric',DTYPE,0,NK,NK,1,0,A(NNZ),2,
     *                     WORK(N+L),LAMBND,WORK,COND,SCALE,DETMAN,
     *                     DETEXP,DIST,SEED)
            END IF
            DO 120 J = 1, NK
               I = I + 1
               L = 2*J + NNZ - 2
               IROW(L) = I
               ICOL(L) = I
               IF (J.LT.NK) THEN
                  IROW(L+1) = I + 1
                  ICOL(L+1) = I
               END IF
  120       CONTINUE
            NNZ = NNZ + 2*NK - 1
  140    CONTINUE
         IF (DTYPE.GE.1) CALL DCOPY(N,WORK(N+1),1,LAMBDA,1)
         NNZ = NNZ - 1
         X = ONE/X02AMF()
         COND = -ONE
         DO 160 I = 1, N
            X = MIN(X,ABS(LAMBDA(I)))
            COND = MAX(COND,ABS(LAMBDA(I)))
  160    CONTINUE
         COND = COND/X
      END IF
*
*     Generate the random diagonal permutations
*
      DO 180 I = 1, N
         IWORK(N+I) = I
  180 CONTINUE
*
      DO 220 I = N + 1, 2*N - 1
  200    CONTINUE
         J = INT((N-I)*Y90TBF(1,SEED)) + I + 1
         IF (J.EQ.I) GO TO 200
         K = IWORK(I)
         IWORK(I) = IWORK(J)
         IWORK(J) = K
  220 CONTINUE
*
      DO 240 I = 1, N
         IWORK(IWORK(N+I)) = I
  240 CONTINUE
*
*     Apply random diagonal permutations
*
      DO 260 K = 1, 2*N - 1
         IROW(K) = IWORK(IROW(K))
         ICOL(K) = IWORK(ICOL(K))
         IF (IROW(K).LT.ICOL(K)) THEN
            I = IROW(K)
            IROW(K) = ICOL(K)
            ICOL(K) = I
         END IF
  260 CONTINUE
*-----------------------------------------------------------------------
*
*     Generate the sparse matrix
*
*-----------------------------------------------------------------------
      CALL Y90RUZ(N,DENS,NNZ,A,IROW,ICOL,IDIMA,SEED,WORK,IWORK)
*
      DENS = HUND*DBLE(2*NNZ-N)/DBLE(N*N)
*
*     End of subroutine Y90RUF
*
      RETURN
      END
