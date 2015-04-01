      SUBROUTINE F08ASX(DIRECT,STOREV,N,K,V,LDV,TAU,T,LDT)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZLARFT(DIRECT,STOREV,N,K,V,LDV,TAU,T,LDT)
C
C  Purpose
C  =======
C
C  ZLARFT forms the triangular factor T of a complex block reflector H
C  of order n, which is defined as a product of k elementary reflectors.
C
C  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
C
C  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
C
C  If STOREV = 'C', the vector which defines the elementary reflector
C  H(i) is stored in the i-th column of the array V, and
C
C     H  =  I - V * T * V'
C
C  If STOREV = 'R', the vector which defines the elementary reflector
C  H(i) is stored in the i-th row of the array V, and
C
C     H  =  I - V' * T * V
C
C  Arguments
C  =========
C
C  DIRECT  (input) CHARACTER*1
C          Specifies the order in which the elementary reflectors are
C          multiplied to form the block reflector:
C          = 'F': H = H(1) H(2) . . . H(k) (Forward)
C          = 'B': H = H(k) . . . H(2) H(1) (Backward)
C
C  STOREV  (input) CHARACTER*1
C          Specifies how the vectors which define the elementary
C          reflectors are stored (see also Further Details):
C          = 'C': columnwise
C          = 'R': rowwise
C
C  N       (input) INTEGER
C          The order of the block reflector H. N >= 0.
C
C  K       (input) INTEGER
C          The order of the triangular factor T (= the number of
C          elementary reflectors). K >= 1.
C
C  V       (input/output) COMPLEX*16 array, dimension
C                               (LDV,K) if STOREV = 'C'
C                               (LDV,N) if STOREV = 'R'
C          The matrix V. See further details.
C
C  LDV     (input) INTEGER
C          The leading dimension of the array V.
C          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
C
C  TAU     (input) COMPLEX*16 array, dimension (K)
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i).
C
C  T       (output) COMPLEX*16 array, dimension (LDT,K)
C          The k by k triangular factor T of the block reflector.
C          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
C          lower triangular. The rest of the array is not used.
C
C  LDT     (input) INTEGER
C          The leading dimension of the array T. LDT >= K.
C
C  Further Details
C  ===============
C
C  The shape of the matrix V and the storage of the vectors which define
C  the H(i) is best illustrated by the following example with n = 5 and
C  k = 3. The elements equal to 1 are not stored; the corresponding
C  array elements are modified but restored on exit. The rest of the
C  array is not used.
C
C  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
C
C               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
C                   ( v1  1    )                     (     1 v2 v2 v2 )
C                   ( v1 v2  1 )                     (        1 v3 v3 )
C                   ( v1 v2 v3 )
C                   ( v1 v2 v3 )
C
C  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
C
C               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
C                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
C                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
C                   (     1 v3 )
C                   (        1 )
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           K, LDT, LDV, N
      CHARACTER         DIRECT, STOREV
C     .. Array Arguments ..
      COMPLEX*16        T(LDT,*), TAU(*), V(LDV,*)
C     .. Local Scalars ..
      COMPLEX*16        VII
      INTEGER           I, J
C     .. External Subroutines ..
      EXTERNAL          F07FRY, ZGEMV, ZTRMV
C     .. Executable Statements ..
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
      IF ((DIRECT.EQ.'F' .OR. DIRECT.EQ.'f')) THEN
         DO 40 I = 1, K
            IF (TAU(I).EQ.ZERO) THEN
C
C              H(i)  =  I
C
               DO 20 J = 1, I
                  T(J,I) = ZERO
   20          CONTINUE
            ELSE
C
C              general case
C
               VII = V(I,I)
               V(I,I) = ONE
               IF ((STOREV.EQ.'C' .OR. STOREV.EQ.'c')) THEN
C
C                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i)
C
                  CALL ZGEMV('Conjugate transpose',N-I+1,I-1,-TAU(I),
     *                       V(I,1),LDV,V(I,I),1,ZERO,T(1,I),1)
               ELSE
C
C                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)'
C
                  IF (I.LT.N) CALL F07FRY(N-I,V(I,I+1),LDV)
                  CALL ZGEMV('No transpose',I-1,N-I+1,-TAU(I),V(1,I),
     *                       LDV,V(I,I),LDV,ZERO,T(1,I),1)
                  IF (I.LT.N) CALL F07FRY(N-I,V(I,I+1),LDV)
               END IF
               V(I,I) = VII
C
C              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
C
               CALL ZTRMV('Upper','No transpose','Non-unit',I-1,T,LDT,
     *                    T(1,I),1)
               T(I,I) = TAU(I)
            END IF
   40    CONTINUE
      ELSE
         DO 80 I = K, 1, -1
            IF (TAU(I).EQ.ZERO) THEN
C
C              H(i)  =  I
C
               DO 60 J = I, K
                  T(J,I) = ZERO
   60          CONTINUE
            ELSE
C
C              general case
C
               IF (I.LT.K) THEN
                  IF ((STOREV.EQ.'C' .OR. STOREV.EQ.'c')) THEN
                     VII = V(N-K+I,I)
                     V(N-K+I,I) = ONE
C
C                    T(i+1:k,i) :=
C                            - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i)
C
                     CALL ZGEMV('Conjugate transpose',N-K+I,K-I,-TAU(I),
     *                          V(1,I+1),LDV,V(1,I),1,ZERO,T(I+1,I),1)
                     V(N-K+I,I) = VII
                  ELSE
                     VII = V(I,N-K+I)
                     V(I,N-K+I) = ONE
C
C                    T(i+1:k,i) :=
C                            - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)'
C
                     CALL F07FRY(N-K+I-1,V(I,1),LDV)
                     CALL ZGEMV('No transpose',K-I,N-K+I,-TAU(I),
     *                          V(I+1,1),LDV,V(I,1),LDV,ZERO,T(I+1,I),1)
                     CALL F07FRY(N-K+I-1,V(I,1),LDV)
                     V(I,N-K+I) = VII
                  END IF
C
C                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
C
                  CALL ZTRMV('Lower','No transpose','Non-unit',K-I,
     *                       T(I+1,I+1),LDT,T(I+1,I),1)
               END IF
               T(I,I) = TAU(I)
            END IF
   80    CONTINUE
      END IF
      RETURN
C
C     End of F08ASX (ZLARFT)
C
      END
