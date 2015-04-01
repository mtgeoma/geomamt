      SUBROUTINE F08HEF(VECT,UPLO,N,KD,AB,LDAB,D,E,Q,LDQ,WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             DSBTRD(VECT,UPLO,N,KD,AB,LDAB,D,E,Q,LDQ,WORK,
     *                  INFO)
C
C  Purpose
C  =======
C
C  DSBTRD reduces a real symmetric band matrix A to symmetric
C  tridiagonal form T by an orthogonal similarity transformation:
C  Q' * A * Q = T.
C
C  Arguments
C  =========
C
C  VECT    (input) CHARACTER*1
C          Specifies whether or not the matrix Q is to be formed.
C          = 'N': do not form Q
C          = 'V': form Q
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
C  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
C          On entry, the upper or lower triangle of the symmetric band
C          matrix A, stored in the first KD+1 rows of the array.  The
C          j-th column of A is stored in the j-th column of the array AB
C          as follows:
C          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
C          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
C          On exit, the diagonal elements of A are overwritten by the
C          diagonal elements of the tridiagonal matrix T; if KD > 0, the
C          elements on the first superdiagonal (if UPLO = 'U') or the
C          first subdiagonal (if UPLO = 'L') are overwritten by the
C          offdiagonal elements of T; the rest of A is overwritten by
C          values generated during the reduction.
C
C  LDAB    (input) INTEGER
C          The leading dimension of the array AB.  LDAB >= KD+1.
C
C  D       (output) DOUBLE PRECISION array, dimension (N)
C          The diagonal elements of the tridiagonal matrix T.
C
C  E       (output) DOUBLE PRECISION array, dimension (N-1)
C          The off-diagonal elements of the tridiagonal matrix T:
C          E(i) = T(i,i+1) if UPLO = 'U'; E(i) = T(i+1,i) if UPLO = 'L'.
C
C  Q       (output) DOUBLE PRECISION array, dimension (LDQ,N)
C          If VECT = 'V', the n by n orthogonal matrix Q.
C          If VECT = 'N', the array Q is not referenced.
C
C  LDQ     (input) INTEGER
C          The leading dimension of the array Q.
C          LDQ >= max(1,N) if VECT = 'V'.
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
      INTEGER           INFO, KD, LDAB, LDQ, N
      CHARACTER         UPLO, VECT
C     .. Array Arguments ..
      DOUBLE PRECISION  AB(LDAB,*), D(*), E(*), Q(LDQ,*), WORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP
      INTEGER           I, J, J1, J2, K, KD1, KDN, L, NR, NRT
      LOGICAL           UPPER, WANTQ
C     .. External Subroutines ..
      EXTERNAL          DROT, F06AAZ, F06QHF, F08HEW, F08HEX, F08HEY,
     *                  F08HEZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     Test the input parameters
C
      WANTQ = (VECT.EQ.'V' .OR. VECT.EQ.'v')
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
      KD1 = KD + 1
      INFO = 0
      IF ( .NOT. WANTQ .AND. .NOT. (VECT.EQ.'N' .OR. VECT.EQ.'n')) THEN
         INFO = -1
      ELSE IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l'))
     *         THEN
         INFO = -2
      ELSE IF (N.LT.0) THEN
         INFO = -3
      ELSE IF (KD.LT.0) THEN
         INFO = -4
      ELSE IF (LDAB.LT.KD1) THEN
         INFO = -6
      ELSE IF ((LDQ.LT.1) .OR. (LDQ.LT.MAX(1,N) .AND. WANTQ)) THEN
         INFO = -10
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08HEF/DSBTRD',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
C     Initialize Q to the unit matrix, if needed
C
      IF (WANTQ) CALL F06QHF('General',N,N,ZERO,ONE,Q,LDQ)
C
C     Wherever possible, plane rotations are generated and applied in
C     vector operations of length NR over the index set J1:J2:KD1.
C
C     The cosines and sines of the plane rotations are stored in the
C     arrays D and WORK.
C
      KDN = MIN(N-1,KD)
      IF (UPPER) THEN
C
         IF (KD.GT.1) THEN
C
C           Reduce to tridiagonal form, working with upper triangle
C
            NR = 0
            J1 = KDN + 2
            J2 = 1
C
            DO 120 I = 1, N - 2
C
C              Reduce i-th row of matrix to tridiagonal form
C
               DO 100 K = KDN + 1, 2, -1
                  J1 = J1 + KDN
                  J2 = J2 + KDN
C
                  IF (NR.GT.0) THEN
C
C                    generate plane rotations to annihilate nonzero
C                    elements which have been created outside the band
C
                     CALL F08HEZ(NR,AB(1,J1-1),KD1*LDAB,WORK(J1),KD1,
     *                           D(J1),KD1)
C
C                    apply rotations from the right
C
                     DO 20 L = 1, KD - 1
                        CALL F08HEY(NR,AB(L+1,J1-1),KD1*LDAB,AB(L,J1),
     *                              KD1*LDAB,D(J1),WORK(J1),KD1)
   20                CONTINUE
                  END IF
C
                  IF (K.GT.2) THEN
                     IF (K.LE.N-I+1) THEN
C
C                       generate plane rotation to annihilate a(i,i+k-1)
C                       within the band
C
                        CALL F08HEW(AB(KD-K+3,I+K-2),AB(KD-K+2,I+K-1),
     *                              D(I+K-1),WORK(I+K-1),TEMP)
                        AB(KD-K+3,I+K-2) = TEMP
C
C                       apply rotation from the right
C
                        CALL DROT(K-3,AB(KD-K+4,I+K-2),1,
     *                            AB(KD-K+3,I+K-1),1,D(I+K-1),
     *                            WORK(I+K-1))
                     END IF
                     NR = NR + 1
                     J1 = J1 - KDN - 1
                  END IF
C
C                 apply plane rotations from both sides to diagonal
C                 blocks
C
                  IF (NR.GT.0) CALL F08HEX(NR,AB(KD1,J1-1),AB(KD1,J1),
     *                                     AB(KD,J1),KD1*LDAB,D(J1),
     *                                     WORK(J1),KD1)
C
C                 apply plane rotations from the left
C
                  DO 40 L = 1, KD - 1
                     IF (J2+L.GT.N) THEN
                        NRT = NR - 1
                     ELSE
                        NRT = NR
                     END IF
                     IF (NRT.GT.0) CALL F08HEY(NRT,AB(KD-L,J1+L),
     *                                  KD1*LDAB,AB(KD-L+1,J1+L),
     *                                  KD1*LDAB,D(J1),WORK(J1),KD1)
   40             CONTINUE
C
                  IF (WANTQ) THEN
C
C                    accumulate product of plane rotations in Q
C
                     DO 60 J = J1, J2, KD1
                        CALL DROT(N,Q(1,J-1),1,Q(1,J),1,D(J),WORK(J))
   60                CONTINUE
                  END IF
C
                  IF (J2+KDN.GT.N) THEN
C
C                    adjust J2 to keep within the bounds of the matrix
C
                     NR = NR - 1
                     J2 = J2 - KDN - 1
                  END IF
C
                  DO 80 J = J1, J2, KD1
C
C                    create nonzero element a(j-1,j+kd) outside the band
C                    and store it in WORK
C
                     WORK(J+KD) = WORK(J)*AB(1,J+KD)
                     AB(1,J+KD) = D(J)*AB(1,J+KD)
   80             CONTINUE
  100          CONTINUE
  120       CONTINUE
         END IF
C
         IF (KD.GT.0) THEN
C
C           copy off-diagonal elements to E
C
            DO 140 I = 1, N - 1
               E(I) = AB(KD,I+1)
  140       CONTINUE
         ELSE
C
C           set E to zero if original matrix was diagonal
C
            DO 160 I = 1, N - 1
               E(I) = ZERO
  160       CONTINUE
         END IF
C
C        copy diagonal elements to D
C
         DO 180 I = 1, N
            D(I) = AB(KD1,I)
  180    CONTINUE
C
      ELSE
C
         IF (KD.GT.1) THEN
C
C           Reduce to tridiagonal form, working with lower triangle
C
            NR = 0
            J1 = KDN + 2
            J2 = 1
C
            DO 300 I = 1, N - 2
C
C              Reduce i-th column of matrix to tridiagonal form
C
               DO 280 K = KDN + 1, 2, -1
                  J1 = J1 + KDN
                  J2 = J2 + KDN
C
                  IF (NR.GT.0) THEN
C
C                    generate plane rotations to annihilate nonzero
C                    elements which have been created outside the band
C
                     CALL F08HEZ(NR,AB(KD1,J1-KD1),KD1*LDAB,WORK(J1),
     *                           KD1,D(J1),KD1)
C
C                    apply plane rotations from one side
C
                     DO 200 L = 1, KD - 1
                        CALL F08HEY(NR,AB(KD1-L,J1-KD1+L),KD1*LDAB,
     *                              AB(KD1-L+1,J1-KD1+L),KD1*LDAB,D(J1),
     *                              WORK(J1),KD1)
  200                CONTINUE
                  END IF
C
                  IF (K.GT.2) THEN
                     IF (K.LE.N-I+1) THEN
C
C                       generate plane rotation to annihilate a(i+k-1,i)
C                       within the band
C
                        CALL F08HEW(AB(K-1,I),AB(K,I),D(I+K-1),
     *                              WORK(I+K-1),TEMP)
                        AB(K-1,I) = TEMP
C
C                       apply rotation from the left
C
                        CALL DROT(K-3,AB(K-2,I+1),LDAB-1,AB(K-1,I+1),
     *                            LDAB-1,D(I+K-1),WORK(I+K-1))
                     END IF
                     NR = NR + 1
                     J1 = J1 - KDN - 1
                  END IF
C
C                 apply plane rotations from both sides to diagonal
C                 blocks
C
                  IF (NR.GT.0) CALL F08HEX(NR,AB(1,J1-1),AB(1,J1),
     *                                     AB(2,J1-1),KD1*LDAB,D(J1),
     *                                     WORK(J1),KD1)
C
C                 apply plane rotations from the right
C
                  DO 220 L = 1, KD - 1
                     IF (J2+L.GT.N) THEN
                        NRT = NR - 1
                     ELSE
                        NRT = NR
                     END IF
                     IF (NRT.GT.0) CALL F08HEY(NRT,AB(L+2,J1-1),
     *                                  KD1*LDAB,AB(L+1,J1),KD1*LDAB,
     *                                  D(J1),WORK(J1),KD1)
  220             CONTINUE
C
                  IF (WANTQ) THEN
C
C                    accumulate product of plane rotations in Q
C
                     DO 240 J = J1, J2, KD1
                        CALL DROT(N,Q(1,J-1),1,Q(1,J),1,D(J),WORK(J))
  240                CONTINUE
                  END IF
C
                  IF (J2+KDN.GT.N) THEN
C
C                    adjust J2 to keep within the bounds of the matrix
C
                     NR = NR - 1
                     J2 = J2 - KDN - 1
                  END IF
C
                  DO 260 J = J1, J2, KD1
C
C                    create nonzero element a(j+kd,j-1) outside the
C                    band and store it in WORK
C
                     WORK(J+KD) = WORK(J)*AB(KD1,J)
                     AB(KD1,J) = D(J)*AB(KD1,J)
  260             CONTINUE
  280          CONTINUE
  300       CONTINUE
         END IF
C
         IF (KD.GT.0) THEN
C
C           copy off-diagonal elements to E
C
            DO 320 I = 1, N - 1
               E(I) = AB(2,I)
  320       CONTINUE
         ELSE
C
C           set E to zero if original matrix was diagonal
C
            DO 340 I = 1, N - 1
               E(I) = ZERO
  340       CONTINUE
         END IF
C
C        copy diagonal elements to D
C
         DO 360 I = 1, N
            D(I) = AB(1,I)
  360    CONTINUE
      END IF
C
      RETURN
C
C     End of F08HEF (DSBTRD)
C
      END
