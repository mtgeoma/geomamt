      SUBROUTINE F07BRF(M,N,KL,KU,AB,LDAB,IPIV,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             ZGBTRF(M,N,KL,KU,AB,LDAB,IPIV,INFO)
C
C  Purpose
C  =======
C
C  ZGBTRF computes an LU factorization of a complex m-by-n band matrix A
C  using partial pivoting with row interchanges.
C
C  This is the blocked version of the algorithm, calling Level 3 BLAS.
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
C  AB      (input/output) COMPLEX array, dimension (LDAB,N)
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
C  elements of U because of fill-in resulting from the row interchanges.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
      INTEGER           NBMAX, LDWORK
      PARAMETER         (NBMAX=64,LDWORK=NBMAX+1)
C     .. Scalar Arguments ..
      INTEGER           INFO, KL, KU, LDAB, M, N
C     .. Array Arguments ..
      COMPLEX*16        AB(LDAB,*)
      INTEGER           IPIV(*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP
      INTEGER           I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP,
     *                  JU, K2, KM, KV, NB, NW
C     .. Local Arrays ..
      COMPLEX*16        WORK13(LDWORK,NBMAX), WORK31(LDWORK,NBMAX)
C     .. External Functions ..
      INTEGER           IZAMAX
      EXTERNAL          IZAMAX
C     .. External Subroutines ..
      EXTERNAL          ZCOPY, ZGEMM, ZGERU, ZSCAL, ZSWAP, ZTRSM,
     *                  F06AAZ, F07ARY, F07BRZ, F07ZAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     KV is the number of superdiagonals in the factor U, allowing for
C     fill-in
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
         CALL F06AAZ('F07BRF/ZGBTRF',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
C
C     Determine the block size for this environment
C
      CALL F07ZAZ(1,'F07BRF',NB,0)
C
C     The block size must not exceed the limit set by the size of the
C     local arrays WORK13 and WORK31.
C
      NB = MIN(NB,NBMAX)
C
      IF (NB.LE.1 .OR. NB.GT.KL) THEN
C
C        Use unblocked code
C
         CALL F07BRZ(M,N,KL,KU,AB,LDAB,IPIV,INFO)
      ELSE
C
C        Use blocked code
C
C        Zero the superdiagonal elements of the work array WORK13
C
         DO 40 J = 1, NB
            DO 20 I = 1, J - 1
               WORK13(I,J) = ZERO
   20       CONTINUE
   40    CONTINUE
C
C        Zero the subdiagonal elements of the work array WORK31
C
         DO 80 J = 1, NB
            DO 60 I = J + 1, NB
               WORK31(I,J) = ZERO
   60       CONTINUE
   80    CONTINUE
C
C        Gaussian elimination with partial pivoting
C
C        Set fill-in elements in columns KU+2 to KV to zero
C
         DO 120 J = KU + 2, MIN(KV,N)
            DO 100 I = KV - J + 2, KL
               AB(I,J) = ZERO
  100       CONTINUE
  120    CONTINUE
C
C        JU is the index of the last column affected by the current
C        stage of the factorization
C
         JU = 1
C
         DO 360 J = 1, MIN(M,N), NB
            JB = MIN(NB,MIN(M,N)-J+1)
C
C           The active part of the matrix is partitioned
C
C              A11   A12   A13
C              A21   A22   A23
C              A31   A32   A33
C
C           Here A11, A21 and A31 denote the current block of JB columns
C           which is about to be factorized. The number of rows in the
C           partitioning are JB, I2, I3 respectively, and the numbers
C           of columns are JB, J2, J3. The superdiagonal elements of A13
C           and the subdiagonal elements of A31 lie outside the band.
C
            I2 = MIN(KL-JB,M-J-JB+1)
            I3 = MIN(JB,M-J-KL+1)
C
C           J2 and J3 are computed after JU has been updated.
C
C           Factorize the current block of JB columns
C
            DO 160 JJ = J, J + JB - 1
C
C              Set fill-in elements in column JJ+KV to zero
C
               IF (JJ+KV.LE.N) THEN
                  DO 140 I = 1, KL
                     AB(I,JJ+KV) = ZERO
  140             CONTINUE
               END IF
C
C              Find pivot and test for singularity. KM is the number of
C              subdiagonal elements in the current column.
C
               KM = MIN(KL,M-JJ)
               JP = IZAMAX(KM+1,AB(KV+1,JJ),1)
               IPIV(JJ) = JP + JJ - J
               IF (AB(KV+JP,JJ).NE.ZERO) THEN
                  JU = MAX(JU,MIN(JJ+KU+JP-1,N))
                  IF (JP.NE.1) THEN
C
C                    Apply interchange to columns J to J+JB-1
C
                     IF (JP+JJ-1.LT.J+KL) THEN
C
                        CALL ZSWAP(JB,AB(KV+1+JJ-J,J),LDAB-1,
     *                             AB(KV+JP+JJ-J,J),LDAB-1)
                     ELSE
C
C                       The interchange affects columns J to JJ-1 of A31
C                       which are stored in the work array WORK31
C
                        CALL ZSWAP(JJ-J,AB(KV+1+JJ-J,J),LDAB-1,
     *                             WORK31(JP+JJ-J-KL,1),LDWORK)
                        CALL ZSWAP(J+JB-JJ,AB(KV+1,JJ),LDAB-1,
     *                             AB(KV+JP,JJ),LDAB-1)
                     END IF
                  END IF
C
C                 Compute multipliers
C
                  CALL ZSCAL(KM,ONE/AB(KV+1,JJ),AB(KV+2,JJ),1)
C
C                 Update trailing submatrix within the band and within
C                 the current block. JM is the index of the last column
C                 which needs to be updated.
C
                  JM = MIN(JU,J+JB-1)
                  IF (JM.GT.JJ) CALL ZGERU(KM,JM-JJ,-ONE,AB(KV+2,JJ),1,
     *                                     AB(KV,JJ+1),LDAB-1,
     *                                     AB(KV+1,JJ+1),LDAB-1)
               ELSE
C
C                 If pivot is zero, set INFO to the index of the pivot
C                 unless a zero pivot has already been found.
C
                  IF (INFO.EQ.0) INFO = JJ
               END IF
C
C              Copy current column of A31 into the work array WORK31
C
               NW = MIN(JJ-J+1,I3)
               IF (NW.GT.0) CALL ZCOPY(NW,AB(KV+KL+1-JJ+J,JJ),1,
     *                                 WORK31(1,JJ-J+1),1)
  160       CONTINUE
            IF (J+JB.LE.N) THEN
C
C              Apply the row interchanges to the other blocks.
C
               J2 = MIN(JU-J+1,KV) - JB
               J3 = MAX(0,JU-J-KV+1)
C
C              Use F07ARY to apply the row interchanges to A12, A22, and
C              A32.
C
               CALL F07ARY(J2,AB(KV+1-JB,J+JB),LDAB-1,1,JB,IPIV(J),1)
C
C              Adjust the pivot indices.
C
               DO 180 I = J, J + JB - 1
                  IPIV(I) = IPIV(I) + J - 1
  180          CONTINUE
C
C              Apply the row interchanges to A13, A23, and A33
C              columnwise.
C
               K2 = J - 1 + JB + J2
               DO 220 I = 1, J3
                  JJ = K2 + I
                  DO 200 II = J + I - 1, J + JB - 1
                     IP = IPIV(II)
                     IF (IP.NE.II) THEN
                        TEMP = AB(KV+1+II-JJ,JJ)
                        AB(KV+1+II-JJ,JJ) = AB(KV+1+IP-JJ,JJ)
                        AB(KV+1+IP-JJ,JJ) = TEMP
                     END IF
  200             CONTINUE
  220          CONTINUE
C
C              Update the relevant part of the trailing submatrix
C
               IF (J2.GT.0) THEN
C
C                 Update A12
C
                  CALL ZTRSM('Left','Lower','No transpose','Unit',JB,J2,
     *                       ONE,AB(KV+1,J),LDAB-1,AB(KV+1-JB,J+JB),
     *                       LDAB-1)
C
                  IF (I2.GT.0) THEN
C
C                    Update A22
C
                     CALL ZGEMM('No transpose','No transpose',I2,J2,JB,
     *                          -ONE,AB(KV+1+JB,J),LDAB-1,
     *                          AB(KV+1-JB,J+JB),LDAB-1,ONE,
     *                          AB(KV+1,J+JB),LDAB-1)
                  END IF
C
                  IF (I3.GT.0) THEN
C
C                    Update A32
C
                     CALL ZGEMM('No transpose','No transpose',I3,J2,JB,
     *                          -ONE,WORK31,LDWORK,AB(KV+1-JB,J+JB),
     *                          LDAB-1,ONE,AB(KV+KL+1-JB,J+JB),LDAB-1)
                  END IF
               END IF
C
               IF (J3.GT.0) THEN
C
C                 Copy the lower triangle of A13 into the work array
C                 WORK13
C
                  DO 260 JJ = 1, J3
                     DO 240 II = JJ, JB
                        WORK13(II,JJ) = AB(II-JJ+1,JJ+J+KV-1)
  240                CONTINUE
  260             CONTINUE
C
C                 Update A13 in the work array
C
                  CALL ZTRSM('Left','Lower','No transpose','Unit',JB,J3,
     *                       ONE,AB(KV+1,J),LDAB-1,WORK13,LDWORK)
C
                  IF (I2.GT.0) THEN
C
C                    Update A23
C
                     CALL ZGEMM('No transpose','No transpose',I2,J3,JB,
     *                          -ONE,AB(KV+1+JB,J),LDAB-1,WORK13,LDWORK,
     *                          ONE,AB(1+JB,J+KV),LDAB-1)
                  END IF
C
                  IF (I3.GT.0) THEN
C
C                    Update A33
C
                     CALL ZGEMM('No transpose','No transpose',I3,J3,JB,
     *                          -ONE,WORK31,LDWORK,WORK13,LDWORK,ONE,
     *                          AB(1+KL,J+KV),LDAB-1)
                  END IF
C
C                 Copy the lower triangle of A13 back into place
C
                  DO 300 JJ = 1, J3
                     DO 280 II = JJ, JB
                        AB(II-JJ+1,JJ+J+KV-1) = WORK13(II,JJ)
  280                CONTINUE
  300             CONTINUE
               END IF
            ELSE
C
C              Adjust the pivot indices.
C
               DO 320 I = J, J + JB - 1
                  IPIV(I) = IPIV(I) + J - 1
  320          CONTINUE
            END IF
C
C           Partially undo the interchanges in the current block to
C           restore the upper triangular form of A31 and copy the upper
C           triangle of A31 back into place
C
            DO 340 JJ = J + JB - 1, J, -1
               JP = IPIV(JJ) - JJ + 1
               IF (JP.NE.1) THEN
C
C                 Apply interchange to columns J to JJ-1
C
                  IF (JP+JJ-1.LT.J+KL) THEN
C
C                    The interchange does not affect A31
C
                     CALL ZSWAP(JJ-J,AB(KV+1+JJ-J,J),LDAB-1,
     *                          AB(KV+JP+JJ-J,J),LDAB-1)
                  ELSE
C
C                    The interchange does affect A31
C
                     CALL ZSWAP(JJ-J,AB(KV+1+JJ-J,J),LDAB-1,
     *                          WORK31(JP+JJ-J-KL,1),LDWORK)
                  END IF
               END IF
C
C              Copy the current column of A31 back into place
C
               NW = MIN(I3,JJ-J+1)
               IF (NW.GT.0) CALL ZCOPY(NW,WORK31(1,JJ-J+1),1,
     *                                 AB(KV+KL+1-JJ+J,JJ),1)
  340       CONTINUE
  360    CONTINUE
      END IF
C
      RETURN
C
C     End of F07BRF (ZGBTRF)
C
      END
