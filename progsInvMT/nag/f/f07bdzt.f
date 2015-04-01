      SUBROUTINE F07BDZ(M,N,KL,KU,AB,LDAB,IPIV,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             DGBTF2(M,N,KL,KU,AB,LDAB,IPIV,INFO)
C
C  Purpose
C  =======
C
C  DGBTF2 computes an LU factorization of a real m-by-n band matrix A
C  using partial pivoting with row interchanges.
C
C  This is the unblocked version of the algorithm, calling Level 2 BLAS.
C
C  Arguments
C  =========
C
C  M       (input) INTEGER
C          The number of rows of the matrix A.  M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix A.  N >= 0.
C
C  KL      (input) INTEGER
C          The number of subdiagonals within the band of A.  KL >= 0.
C
C  KU      (input) INTEGER
C          The number of superdiagonals within the band of A.  KU >= 0.
C
C  AB      (input/output) REAL array, dimension (LDAB,N)
C          On entry, the matrix A in band storage, in rows KL+1 to
C          2*KL+KU+1; rows 1 to KL of the array need not be set.
C          The j-th column of A is stored in the j-th column of the
C          array AB as follows:
C          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
C
C          On exit, details of the factorization: U is stored as an
C          upper triangular band matrix with KL+KU superdiagonals in
C          rows 1 to KL+KU+1, and the multipliers used during the
C          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
C          See below for further details.
C
C  LDAB    (input) INTEGER
C          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
C
C  IPIV    (output) INTEGER array, dimension (min(M,N))
C          The pivot indices; for 1 <= i <= min(M,N), row i of the
C          matrix was interchanged with row IPIV(i).
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument had an illegal value
C          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
C               has been completed, but the factor U is exactly
C               singular, and division by zero will occur if it is used
C               to solve a system of equations.
C
C  Further Details
C  ===============
C
C  The band storage scheme is illustrated by the following example, when
C  M = N = 6, KL = 2, KU = 1:
C
C  On entry:                       On exit:
C
C      *    *    *    +    +    +       *    *    *   u14  u25  u36
C      *    *    +    +    +    +       *    *   u13  u24  u35  u46
C      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
C     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
C     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
C     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
C
C  Array elements marked * are not used by the routine; elements marked
C  + need not be set on entry, but are required by the routine to store
C  elements of U, because of fill-in resulting from the row
C  interchanges.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, KL, KU, LDAB, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  AB(LDAB,*)
      INTEGER           IPIV(*)
C     .. Local Scalars ..
      INTEGER           I, J, JP, JU, KM, KV
C     .. External Functions ..
      INTEGER           IDAMAX
      EXTERNAL          IDAMAX
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, DGER, DSCAL, DSWAP
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     KV is the number of superdiagonals in the factor U, allowing for
C     fill-in.
C
      KV = KU + KL
C
C     Test the input parameters.
C
      INFO = 0
      IF (M.LT.0) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (KL.LT.0) THEN
         INFO = -3
      ELSE IF (KU.LT.0) THEN
         INFO = -4
      ELSE IF (LDAB.LT.KL+KV+1) THEN
         INFO = -6
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07BDZ/DGBTF2',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
C
C     Gaussian elimination with partial pivoting
C
C     Set fill-in elements in columns KU+2 to KV to zero.
C
      DO 40 J = KU + 2, MIN(KV,N)
         DO 20 I = KV - J + 2, KL
            AB(I,J) = ZERO
   20    CONTINUE
   40 CONTINUE
C
C     JU is the index of the last column affected by the current stage
C     of the factorization.
C
      JU = 1
C
      DO 80 J = 1, MIN(M,N)
C
C        Set fill-in elements in column J+KV to zero.
C
         IF (J+KV.LE.N) THEN
            DO 60 I = 1, KL
               AB(I,J+KV) = ZERO
   60       CONTINUE
         END IF
C
C        Find pivot and test for singularity. KM is the number of
C        subdiagonal elements in the current column.
C
         KM = MIN(KL,M-J)
         JP = IDAMAX(KM+1,AB(KV+1,J),1)
         IPIV(J) = JP + J - 1
         IF (AB(KV+JP,J).NE.ZERO) THEN
            JU = MAX(JU,MIN(J+KU+JP-1,N))
C
C           Apply interchange to columns J to JU.
C
            IF (JP.NE.1) CALL DSWAP(JU-J+1,AB(KV+JP,J),LDAB-1,AB(KV+1,J)
     *                              ,LDAB-1)
C
            IF (KM.GT.0) THEN
C
C              Compute multipliers.
C
               CALL DSCAL(KM,ONE/AB(KV+1,J),AB(KV+2,J),1)
C
C              Update trailing submatrix within the band.
C
               IF (JU.GT.J) CALL DGER(KM,JU-J,-ONE,AB(KV+2,J),1,
     *                                AB(KV,J+1),LDAB-1,AB(KV+1,J+1),
     *                                LDAB-1)
            END IF
         ELSE
C
C           If pivot is zero, set INFO to the index of the pivot
C           unless a zero pivot has already been found.
C
            IF (INFO.EQ.0) INFO = J
         END IF
   80 CONTINUE
      RETURN
C
C     End of F07BDZ (DGBTF2)
C
      END
