      SUBROUTINE F07HDF(UPLO,N,KD,AB,LDAB,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             DPBTRF(UPLO,N,KD,AB,LDAB,INFO)
C
C  Purpose
C  =======
C
C  DPBTRF computes the Cholesky factorization of a real symmetric
C  positive definite band matrix A.
C
C  The factorization has the form
C     A = U' * U ,  if UPLO = 'U', or
C     A = L  * L',  if UPLO = 'L',
C  where U is an upper triangular matrix and L is lower triangular.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the upper or lower triangular part of the
C          symmetric matrix A is stored:
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  KD      (input) INTEGER
C          The number of super-diagonals of the matrix A if UPLO = 'U',
C          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
C
C  AB      (input/output) REAL array, dimension (LDAB,N)
C          On entry, the upper or lower triangle of the symmetric band
C          matrix A, stored in the first KD+1 rows of the array.  The
C          j-th column of A is stored in the j-th column of the array AB
C          as follows:
C          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
C          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
C
C          On exit, if INFO = 0, the triangular factor U or L from the
C          Cholesky factorization A = U'*U or A = L*L' of the band
C          matrix A, in the same storage format as A.
C
C  LDAB    (input) INTEGER
C          The leading dimension of the array AB.  LDAB >= KD+1.
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = k, the leading minor of order k is not
C               positive definite, and the factorization could not be
C               completed.
C
C  Further Details
C  ===============
C
C  The band storage scheme is illustrated by the following example, when
C  N = 6, KD = 2, and UPLO = 'U':
C
C  On entry:                       On exit:
C
C      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
C      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
C     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
C
C  Similarly, if UPLO = 'L' the format of A is as follows:
C
C  On entry:                       On exit:
C
C     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
C     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
C     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
C
C  Array elements marked * are not used by the routine.
C
C  Contributed by
C  Peter Mayes and Giuseppe Radicati, IBM ECSEC, Rome, March 23, 1989
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
      INTEGER           NBMAX, LDWORK
      PARAMETER         (NBMAX=32,LDWORK=NBMAX+1)
C     .. Scalar Arguments ..
      INTEGER           INFO, KD, LDAB, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  AB(LDAB,*)
C     .. Local Scalars ..
      INTEGER           I, I2, I3, IB, II, J, JJ, NB
C     .. Local Arrays ..
      DOUBLE PRECISION  WORK(LDWORK,NBMAX)
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07FDZ, F07HDZ, F07ZAZ, DGEMM, DSYRK,
     *                  DTRSM
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF (( .NOT. (UPLO.EQ.'U' .OR. UPLO.EQ.'u'))
     *    .AND. ( .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l'))) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (KD.LT.0) THEN
         INFO = -3
      ELSE IF (LDAB.LT.KD+1) THEN
         INFO = -5
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07HDF/DPBTRF',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
C     Determine the block size for this environment
C
      CALL F07ZAZ(1,'F07HDF',NB,0)
C
C     The block size must not exceed the semi-bandwidth KD, and must not
C     exceed the limit set by the size of the local array WORK.
C
      NB = MIN(NB,NBMAX)
C
      IF (NB.LE.1 .OR. NB.GT.KD) THEN
C
C        Use unblocked code
C
         CALL F07HDZ(UPLO,N,KD,AB,LDAB,INFO)
      ELSE
C
C        Use blocked code
C
         IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
C
C           Compute the Cholesky factorization of a symmetric band
C           matrix, given the upper triangle of the matrix in band
C           storage.
C
C           Zero the upper triangle of the work array.
C
            DO 40 J = 1, NB
               DO 20 I = 1, J - 1
                  WORK(I,J) = ZERO
   20          CONTINUE
   40       CONTINUE
C
C           Process the band matrix one diagonal block at a time.
C
            DO 140 I = 1, N, NB
               IB = MIN(NB,N-I+1)
C
C              Factorize the diagonal block
C
               CALL F07FDZ(UPLO,IB,AB(KD+1,I),LDAB-1,II)
               IF (II.NE.0) THEN
                  INFO = I + II - 1
                  GO TO 300
               END IF
               IF (I+IB.LE.N) THEN
C
C                 Update the relevant part of the trailing submatrix.
C                 If A11 denotes the diagonal block which has just been
C                 factorized, then we need to update the remaining
C                 blocks in the diagram:
C
C                    A11   A12   A13
C                          A22   A23
C                                A33
C
C                 The numbers of rows and columns in the partitioning
C                 are IB, I2, I3 respectively. The blocks A12, A22 and
C                 A23 are empty if IB = KD. The upper triangle of A13
C                 lies outside the band.
C
                  I2 = MIN(KD-IB,N-I-IB+1)
                  I3 = MIN(IB,N-I-KD+1)
C
                  IF (I2.GT.0) THEN
C
C                    Update A12
C
                     CALL DTRSM('Left','Upper','Transpose','Non-unit',
     *                          IB,I2,ONE,AB(KD+1,I),LDAB-1,
     *                          AB(KD+1-IB,I+IB),LDAB-1)
C
C                    Update A22
C
                     CALL DSYRK('Upper','Transpose',I2,IB,-ONE,
     *                          AB(KD+1-IB,I+IB),LDAB-1,ONE,
     *                          AB(KD+1,I+IB),LDAB-1)
                  END IF
C
                  IF (I3.GT.0) THEN
C
C                    Copy the lower triangle of A13 into the work array.
C
                     DO 80 JJ = 1, I3
                        DO 60 II = JJ, IB
                           WORK(II,JJ) = AB(II-JJ+1,JJ+I+KD-1)
   60                   CONTINUE
   80                CONTINUE
C
C                    Update A13 (in the work array).
C
                     CALL DTRSM('Left','Upper','Transpose','Non-unit',
     *                          IB,I3,ONE,AB(KD+1,I),LDAB-1,WORK,LDWORK)
C
C                    Update A23
C
                     IF (I2.GT.0) CALL DGEMM('Transpose','No Transpose',
     *                                       I2,I3,IB,-ONE,
     *                                       AB(KD+1-IB,I+IB),LDAB-1,
     *                                       WORK,LDWORK,ONE,
     *                                       AB(1+IB,I+KD),LDAB-1)
C
C                    Update A33
C
                     CALL DSYRK('Upper','Transpose',I3,IB,-ONE,WORK,
     *                          LDWORK,ONE,AB(KD+1,I+KD),LDAB-1)
C
C                    Copy the lower triangle of A13 back into place.
C
                     DO 120 JJ = 1, I3
                        DO 100 II = JJ, IB
                           AB(II-JJ+1,JJ+I+KD-1) = WORK(II,JJ)
  100                   CONTINUE
  120                CONTINUE
                  END IF
               END IF
  140       CONTINUE
         ELSE
C
C           Compute the Cholesky factorization of a symmetric band
C           matrix, given the lower triangle of the matrix in band
C           storage.
C
C           Zero the lower triangle of the work array.
C
            DO 180 J = 1, NB
               DO 160 I = J + 1, NB
                  WORK(I,J) = ZERO
  160          CONTINUE
  180       CONTINUE
C
C           Process the band matrix one diagonal block at a time.
C
            DO 280 I = 1, N, NB
               IB = MIN(NB,N-I+1)
C
C              Factorize the diagonal block
C
               CALL F07FDZ(UPLO,IB,AB(1,I),LDAB-1,II)
               IF (II.NE.0) THEN
                  INFO = I + II - 1
                  GO TO 300
               END IF
               IF (I+IB.LE.N) THEN
C
C                 Update the relevant part of the trailing submatrix.
C                 If A11 denotes the diagonal block which has just been
C                 factorized, then we need to update the remaining
C                 blocks in the diagram:
C
C                    A11
C                    A21   A22
C                    A31   A32   A33
C
C                 The numbers of rows and columns in the partitioning
C                 are IB, I2, I3 respectively. The blocks A21, A22 and
C                 A32 are empty if IB = KD. The lower triangle of A31
C                 lies outside the band.
C
                  I2 = MIN(KD-IB,N-I-IB+1)
                  I3 = MIN(IB,N-I-KD+1)
C
                  IF (I2.GT.0) THEN
C
C                    Update A21
C
                     CALL DTRSM('Right','Lower','Transpose','Non-unit',
     *                          I2,IB,ONE,AB(1,I),LDAB-1,AB(1+IB,I),
     *                          LDAB-1)
C
C                    Update A22
C
                     CALL DSYRK('Lower','No Transpose',I2,IB,-ONE,
     *                          AB(1+IB,I),LDAB-1,ONE,AB(1,I+IB),LDAB-1)
                  END IF
C
                  IF (I3.GT.0) THEN
C
C                    Copy the upper triangle of A31 into the work array.
C
                     DO 220 JJ = 1, IB
                        DO 200 II = 1, MIN(JJ,I3)
                           WORK(II,JJ) = AB(KD+1-JJ+II,JJ+I-1)
  200                   CONTINUE
  220                CONTINUE
C
C                    Update A31 (in the work array).
C
                     CALL DTRSM('Right','Lower','Transpose','Non-unit',
     *                          I3,IB,ONE,AB(1,I),LDAB-1,WORK,LDWORK)
C
C                    Update A32
C
                     IF (I2.GT.0) CALL DGEMM('No transpose','Transpose',
     *                                       I3,I2,IB,-ONE,WORK,LDWORK,
     *                                       AB(1+IB,I),LDAB-1,ONE,
     *                                       AB(1+KD-IB,I+IB),LDAB-1)
C
C                    Update A33
C
                     CALL DSYRK('Lower','No Transpose',I3,IB,-ONE,WORK,
     *                          LDWORK,ONE,AB(1,I+KD),LDAB-1)
C
C                    Copy the upper triangle of A31 back into place.
C
                     DO 260 JJ = 1, IB
                        DO 240 II = 1, MIN(JJ,I3)
                           AB(KD+1-JJ+II,JJ+I-1) = WORK(II,JJ)
  240                   CONTINUE
  260                CONTINUE
                  END IF
               END IF
  280       CONTINUE
         END IF
      END IF
      RETURN
C
  300 CONTINUE
      RETURN
C
C     End of F07HDF (DPBTRF)
C
      END
