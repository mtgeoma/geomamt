      SUBROUTINE F08NEY(N,K,NB,A,LDA,TAU,T,LDT,Y,LDY)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DLAHRD(N,K,NB,A,LDA,TAU,T,LDT,Y,LDY)
C
C  Purpose
C  =======
C
C  DLAHRD reduces the first NB columns of a real general n-by-(n-k+1)
C  matrix A so that elements below the k-th subdiagonal are zero. The
C  reduction is performed by an orthogonal similarity transformation
C  Q' * A * Q. The routine returns the matrices V and T which determine
C  Q as a block reflector I - V*T*V', and also the matrix Y = A * V * T.
C
C  This is an auxiliary routine called by DGEHRD.
C
C  Arguments
C  =========
C
C  N       (input) INTEGER
C          The order of the matrix A.
C
C  K       (input) INTEGER
C          The offset for the reduction. Elements below the k-th
C          subdiagonal in the first NB columns are reduced to zero.
C
C  NB      (input) INTEGER
C          The number of columns to be reduced.
C
C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N-K+1)
C          On entry, the n-by-(n-k+1) general matrix A.
C          On exit, the elements on and above the k-th subdiagonal in
C          the first NB columns are overwritten with the corresponding
C          elements of the reduced matrix; the elements below the k-th
C          subdiagonal, with the array TAU, represent the matrix Q as a
C          product of elementary reflectors. The other columns of A are
C          unchanged. See Further Details.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  TAU     (output) DOUBLE PRECISION array, dimension (NB)
C          The scalar factors of the elementary reflectors. See Further
C          Details.
C
C  T       (output) DOUBLE PRECISION array, dimension (NB,NB)
C          The upper triangular matrix T.
C
C  LDT     (input) INTEGER
C          The leading dimension of the array T.  LDT >= NB.
C
C  Y       (output) DOUBLE PRECISION array, dimension (LDY,NB)
C          The n-by-nb matrix Y.
C
C  LDY     (input) INTEGER
C          The leading dimension of the array Y. LDY >= N.
C
C  Further Details
C  ===============
C
C  The matrix Q is represented as a product of nb elementary reflectors
C
C     Q = H(1) H(2) . . . H(nb).
C
C  Each H(i) has the form
C
C     H(i) = I - tau * v * v'
C
C  where tau is a real scalar, and v is a real vector with
C  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in
C  A(i+k+1:n,i), and tau in TAU(i).
C
C  The elements of the vectors v together form the (n-k+1)-by-nb matrix
C  V which is needed, with T and Y, to apply the transformation to the
C  unreduced part of the matrix, using an update of the form:
C  A := (I - V*T*V') * (A - Y*V').
C
C  The contents of A on exit are illustrated by the following example
C  with n = 7, k = 3 and nb = 2:
C
C     ( a   h   a   a   a )
C     ( a   h   a   a   a )
C     ( a   h   a   a   a )
C     ( h   h   a   a   a )
C     ( v1  h   a   a   a )
C     ( v1  v2  a   a   a )
C     ( v1  v2  a   a   a )
C
C  where a denotes an element of the original matrix A, h denotes a
C  modified element of the upper Hessenberg matrix H, and vi denotes an
C  element of the vector defining H(i).
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           K, LDA, LDT, LDY, N, NB
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), T(LDT,NB), TAU(NB), Y(LDY,NB)
C     .. Local Scalars ..
      DOUBLE PRECISION  EI
      INTEGER           I
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DSCAL, DTRMV, F08AEV
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
C     Quick return if possible
C
      IF (N.LE.1) RETURN
C
      DO 20 I = 1, NB
         IF (I.GT.1) THEN
C
C           Update A(1:n,i)
C
C           Compute i-th column of A - Y * V'
C
            CALL DGEMV('No transpose',N,I-1,-ONE,Y,LDY,A(K+I-1,1),LDA,
     *                 ONE,A(1,I),1)
C
C           Apply I - V * T' * V' to this column (call it b) from the
C           left, using the last column of T as workspace
C
C           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
C                    ( V2 )             ( b2 )
C
C           where V1 is unit lower triangular
C
C           w := V1' * b1
C
            CALL DCOPY(I-1,A(K+1,I),1,T(1,NB),1)
            CALL DTRMV('Lower','Transpose','Unit',I-1,A(K+1,1),LDA,
     *                 T(1,NB),1)
C
C           w := w + V2'*b2
C
            CALL DGEMV('Transpose',N-K-I+1,I-1,ONE,A(K+I,1),LDA,A(K+I,I)
     *                 ,1,ONE,T(1,NB),1)
C
C           w := T'*w
C
            CALL DTRMV('Upper','Transpose','Non-unit',I-1,T,LDT,T(1,NB),
     *                 1)
C
C           b2 := b2 - V2*w
C
            CALL DGEMV('No transpose',N-K-I+1,I-1,-ONE,A(K+I,1),LDA,
     *                 T(1,NB),1,ONE,A(K+I,I),1)
C
C           b1 := b1 - V1*w
C
            CALL DTRMV('Lower','No transpose','Unit',I-1,A(K+1,1),LDA,
     *                 T(1,NB),1)
            CALL DAXPY(I-1,-ONE,T(1,NB),1,A(K+1,I),1)
C
            A(K+I-1,I-1) = EI
         END IF
C
C        Generate the elementary reflector H(i) to annihilate
C        A(k+i+1:n,i)
C
         CALL F08AEV(N-K-I+1,A(K+I,I),A(MIN(K+I+1,N),I),1,TAU(I))
         EI = A(K+I,I)
         A(K+I,I) = ONE
C
C        Compute  Y(1:n,i)
C
         CALL DGEMV('No transpose',N,N-K-I+1,ONE,A(1,I+1),LDA,A(K+I,I),
     *              1,ZERO,Y(1,I),1)
         CALL DGEMV('Transpose',N-K-I+1,I-1,ONE,A(K+I,1),LDA,A(K+I,I),1,
     *              ZERO,T(1,I),1)
         CALL DGEMV('No transpose',N,I-1,-ONE,Y,LDY,T(1,I),1,ONE,Y(1,I),
     *              1)
         CALL DSCAL(N,TAU(I),Y(1,I),1)
C
C        Compute T(1:i,i)
C
         CALL DSCAL(I-1,-TAU(I),T(1,I),1)
         CALL DTRMV('Upper','No transpose','Non-unit',I-1,T,LDT,T(1,I),
     *              1)
         T(I,I) = TAU(I)
C
   20 CONTINUE
      A(K+NB,NB) = EI
C
      RETURN
C
C     End of F08NEY (DLAHRD)
C
      END
