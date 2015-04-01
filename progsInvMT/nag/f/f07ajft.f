      SUBROUTINE F07AJF(N,A,LDA,IPIV,WORK,LWORK,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             DGETRI(N,A,LDA,IPIV,WORK,LWORK,INFO)
C
C  Purpose
C  =======
C
C  DGETRI computes the inverse of a matrix using the LU factorization
C  computed by F07ADF.
C
C  This method inverts U and then computes inv(A) by solving the system
C  inv(A)*L = inv(U) for inv(A).
C
C  Arguments
C  =========
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  A       (input/output) REAL array, dimension (LDA,N)
C          On entry, the factors L and U from the factorization
C          A = P*L*U as computed by F07ADF.
C          On exit, if INFO = 0, the inverse of the original matrix A.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  IPIV    (input) INTEGER array, dimension (N)
C          The pivot indices from F07ADF; for 1<=i<=N, row i of the
C          matrix was interchanged with row IPIV(i).
C
C  WORK    (workspace) REAL array, dimension (LWORK)
C          If INFO returns 0, then WORK(1) returns the minimum
C          value of LWORK required for optimal performance.
C
C  LWORK   (input) INTEGER
C          The dimension of the array WORK.  LWORK >= max(N,1).
C          For optimal performance LWORK should be at least N*NB,
C          where NB is the optimal blocksize returned by F07ZAZ.
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = k, U(k,k) is exactly zero; the matrix is
C               singular and its inverse could not be computed.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LWORK, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), WORK(LWORK)
      INTEGER           IPIV(*)
C     .. Local Scalars ..
      INTEGER           I, IWS, J, JB, JJ, JP, LDWORK, NB, NBMIN, NN
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07TJF, F07ZAZ, DGEMM, DGEMV, DSWAP,
     *                  DTRSM
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      WORK(1) = MAX(N,1)
      IF (N.LT.0) THEN
         INFO = -1
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -3
      ELSE IF (LWORK.LT.N) THEN
         INFO = -6
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07AJF/DGETRI',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) THEN
         WORK(1) = 1
         RETURN
      END IF
C
C     Form inv(U).  If INFO > 0 from F07TJF, then U is singular,
C     and the inverse is not computed.
C
      CALL F07TJF('Upper','Non-unit',N,A,LDA,INFO)
      IF (INFO.GT.0) RETURN
C
C     Determine the block size for this environment.
C
      CALL F07ZAZ(1,'F07AJF',NB,0)
      IF (NB.LE.1) NB = N
C
      IF (NB.LT.N) THEN
         LDWORK = N
C
C        Determine if workspace is large enough for blocked code
C
         IWS = N*NB
         IF (LWORK.LT.IWS) THEN
C
C           Not enough workspace has been supplied to use the optimal
C           value of NB: determine the minimum value of NB, and reduce
C           NB or force use of unblocked code
C
            CALL F07ZAZ(2,'F07AJF',NBMIN,0)
            NBMIN = MAX(2,NBMIN)
C
            IF (LWORK.GE.N*NBMIN) THEN
               NB = LWORK/N
            ELSE
               NB = N
            END IF
         END IF
      ELSE
         IWS = N
      END IF
C
C     Solve the equation inv(A)*L = inv(U) for inv(A).
C
      IF (NB.GE.N) THEN
C
C        Use unblocked code.
C
         DO 40 J = N, 1, -1
C
C           Copy current column of L to WORK and replace with zeros.
C
            DO 20 I = J + 1, N
               WORK(I) = A(I,J)
               A(I,J) = ZERO
   20       CONTINUE
C
C           Compute current column of inv(A).
C
            IF (J.LT.N) CALL DGEMV('No transpose',N,N-J,-ONE,A(1,J+1),
     *                             LDA,WORK(J+1),1,ONE,A(1,J),1)
   40    CONTINUE
      ELSE
C
C        Use blocked code.
C
         NN = ((N-1)/NB)*NB + 1
         DO 100 J = NN, 1, -NB
            JB = MIN(NB,N-J+1)
C
C           Copy current block column of L to WORK and replace with
C           zeros.
C
            DO 80 JJ = J, J + JB - 1
               DO 60 I = JJ + 1, N
                  WORK(I+(JJ-J)*LDWORK) = A(I,JJ)
                  A(I,JJ) = ZERO
   60          CONTINUE
   80       CONTINUE
C
C           Compute current block column of inv(A).
C
            IF (J+JB.LE.N) CALL DGEMM('No transpose','No transpose',N,
     *                                JB,N-J-JB+1,-ONE,A(1,J+JB),LDA,
     *                                WORK(J+JB),LDWORK,ONE,A(1,J),LDA)
            CALL DTRSM('Right','Lower','No transpose','Unit',N,JB,ONE,
     *                 WORK(J),LDWORK,A(1,J),LDA)
  100    CONTINUE
      END IF
C
C     Apply column interchanges.
C
      DO 120 J = N - 1, 1, -1
         JP = IPIV(J)
         IF (JP.NE.J) CALL DSWAP(N,A(1,J),1,A(1,JP),1)
  120 CONTINUE
C
      WORK(1) = IWS
      RETURN
C
C     End of F07AJF (DGETRI)
C
      END
