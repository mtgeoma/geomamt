      SUBROUTINE F08FEY(UPLO,N,NB,A,LDA,E,TAU,W,LDW)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DLATRD(UPLO,N,NB,A,LDA,E,TAU,W,LDW)
C
C  Purpose
C  =======
C
C  DLATRD reduces NB rows and columns of a real symmetric matrix A to
C  symmetric tridiagonal form by an orthogonal similarity
C  transformation Q' * A * Q, and returns the matrices V and W which are
C  needed to apply the transformation to the unreduced part of A.
C
C  If UPLO = 'U', DLATRD reduces the last NB rows and columns of a
C  matrix, of which the upper triangle is supplied;
C  if UPLO = 'L', DLATRD reduces the first NB rows and columns of a
C  matrix, of which the lower triangle is supplied.
C
C  This is an auxiliary routine called by DSYTRD.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER
C          Specifies whether the upper or lower triangular part of the
C          symmetric matrix A is stored:
C          = 'U': Upper triangular
C          = 'L': Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.
C
C  NB      (input) INTEGER
C          The number of rows and columns to be reduced.
C
C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
C          n-by-n upper triangular part of A contains the upper
C          triangular part of the matrix A, and the strictly lower
C          triangular part of A is not referenced.  If UPLO = 'L', the
C          leading n-by-n lower triangular part of A contains the lower
C          triangular part of the matrix A, and the strictly upper
C          triangular part of A is not referenced.
C          On exit:
C          if UPLO = 'U', the last NB columns have been reduced to
C            tridiagonal form, with the diagonal elements overwriting
C            the diagonal elements of A; the elements above the diagonal
C            with the array TAU, represent the orthogonal matrix Q as a
C            product of elementary reflectors;
C          if UPLO = 'L', the first NB columns have been reduced to
C            tridiagonal form, with the diagonal elements overwriting
C            the diagonal elements of A; the elements below the diagonal
C            with the array TAU, represent the  orthogonal matrix Q as a
C            product of elementary reflectors.
C          See Further Details.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= (1,N).
C
C  E       (output) DOUBLE PRECISION array, dimension (N-1)
C          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
C          elements of the last NB columns of the reduced matrix;
C          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of
C          the first NB columns of the reduced matrix.
C
C  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
C          The scalar factors of the elementary reflectors, stored in
C          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
C          See Further Details.
C
C  W       (output) DOUBLE PRECISION array, dimension (LDW,NB)
C          The n-by-nb matrix W required to update the unreduced part
C          of A.
C
C  LDW     (input) INTEGER
C          The leading dimension of the array W. LDW >= max(1,N).
C
C  Further Details
C  ===============
C
C  If UPLO = 'U', the matrix Q is represented as a product of elementary
C  reflectors
C
C     Q = H(n) H(n-1) . . . H(n-nb+1).
C
C  Each H(i) has the form
C
C     H(i) = I - tau * v * v'
C
C  where tau is a real scalar, and v is a real vector with
C  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
C  and tau in TAU(i-1).
C
C  If UPLO = 'L', the matrix Q is represented as a product of elementary
C  reflectors
C
C     Q = H(1) H(2) . . . H(nb).
C
C  Each H(i) has the form
C
C     H(i) = I - tau * v * v'
C
C  where tau is a real scalar, and v is a real vector with
C  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
C  and tau in TAU(i).
C
C  The elements of the vectors v together form the n-by-nb matrix V
C  which is needed, with W, to apply the transformation to the unreduced
C  part of the matrix, using a symmetric rank-2k update of the form:
C  A := A - V*W' - W*V'.
C
C  The contents of A on exit are illustrated by the following examples
C  with n = 5 and nb = 2:
C
C  if UPLO = 'U':                       if UPLO = 'L':
C
C    (  a   a   a   v4  v5 )              (  d                  )
C    (      a   a   v4  v5 )              (  1   d              )
C    (          a   1   v5 )              (  v1  1   a          )
C    (              d   1  )              (  v1  v2  a   a      )
C    (                  d  )              (  v1  v2  a   a   a  )
C
C  where d denotes a diagonal element of the reduced matrix, a denotes
C  an element of the original matrix that is unchanged, and vi denotes
C  an element of the vector defining H(i).
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, HALF
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0,HALF=0.5D+0)
C     .. Scalar Arguments ..
      INTEGER           LDA, LDW, N, NB
      CHARACTER         UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), E(*), TAU(*), W(LDW,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA
      INTEGER           I, IW
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DGEMV, DSCAL, DSYMV, F08AEV
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
C     Quick return if possible
C
      IF (N.LE.0) RETURN
C
      IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
C
C        Reduce last NB columns of upper triangle
C
         DO 20 I = N, N - NB + 1, -1
            IW = I - N + NB
            IF (I.LT.N) THEN
C
C              Update A(1:i,i)
C
               CALL DGEMV('No transpose',I,N-I,-ONE,A(1,I+1),LDA,
     *                    W(I,IW+1),LDW,ONE,A(1,I),1)
               CALL DGEMV('No transpose',I,N-I,-ONE,W(1,IW+1),LDW,
     *                    A(I,I+1),LDA,ONE,A(1,I),1)
            END IF
            IF (I.GT.1) THEN
C
C              Generate elementary reflector H(i) to annihilate
C              A(1:i-2,i)
C
               CALL F08AEV(I-1,A(I-1,I),A(1,I),1,TAU(I-1))
               E(I-1) = A(I-1,I)
               A(I-1,I) = ONE
C
C              Compute W(1:i-1,i)
C
               CALL DSYMV('Upper',I-1,ONE,A,LDA,A(1,I),1,ZERO,W(1,IW),1)
               IF (I.LT.N) THEN
                  CALL DGEMV('Transpose',I-1,N-I,ONE,W(1,IW+1),LDW,
     *                       A(1,I),1,ZERO,W(I+1,IW),1)
                  CALL DGEMV('No transpose',I-1,N-I,-ONE,A(1,I+1),LDA,
     *                       W(I+1,IW),1,ONE,W(1,IW),1)
                  CALL DGEMV('Transpose',I-1,N-I,ONE,A(1,I+1),LDA,A(1,I)
     *                       ,1,ZERO,W(I+1,IW),1)
                  CALL DGEMV('No transpose',I-1,N-I,-ONE,W(1,IW+1),LDW,
     *                       W(I+1,IW),1,ONE,W(1,IW),1)
               END IF
               CALL DSCAL(I-1,TAU(I-1),W(1,IW),1)
               ALPHA = -HALF*TAU(I-1)*DDOT(I-1,W(1,IW),1,A(1,I),1)
               CALL DAXPY(I-1,ALPHA,A(1,I),1,W(1,IW),1)
            END IF
C
   20    CONTINUE
      ELSE
C
C        Reduce first NB columns of lower triangle
C
         DO 40 I = 1, NB
C
C           Update A(i:n,i)
C
            CALL DGEMV('No transpose',N-I+1,I-1,-ONE,A(I,1),LDA,W(I,1),
     *                 LDW,ONE,A(I,I),1)
            CALL DGEMV('No transpose',N-I+1,I-1,-ONE,W(I,1),LDW,A(I,1),
     *                 LDA,ONE,A(I,I),1)
            IF (I.LT.N) THEN
C
C              Generate elementary reflector H(i) to annihilate
C              A(i+2:n,i)
C
               CALL F08AEV(N-I,A(I+1,I),A(MIN(I+2,N),I),1,TAU(I))
               E(I) = A(I+1,I)
               A(I+1,I) = ONE
C
C              Compute W(i+1:n,i)
C
               CALL DSYMV('Lower',N-I,ONE,A(I+1,I+1),LDA,A(I+1,I),1,
     *                    ZERO,W(I+1,I),1)
               CALL DGEMV('Transpose',N-I,I-1,ONE,W(I+1,1),LDW,A(I+1,I),
     *                    1,ZERO,W(1,I),1)
               CALL DGEMV('No transpose',N-I,I-1,-ONE,A(I+1,1),LDA,
     *                    W(1,I),1,ONE,W(I+1,I),1)
               CALL DGEMV('Transpose',N-I,I-1,ONE,A(I+1,1),LDA,A(I+1,I),
     *                    1,ZERO,W(1,I),1)
               CALL DGEMV('No transpose',N-I,I-1,-ONE,W(I+1,1),LDW,
     *                    W(1,I),1,ONE,W(I+1,I),1)
               CALL DSCAL(N-I,TAU(I),W(I+1,I),1)
               ALPHA = -HALF*TAU(I)*DDOT(N-I,W(I+1,I),1,A(I+1,I),1)
               CALL DAXPY(N-I,ALPHA,A(I+1,I),1,W(I+1,I),1)
            END IF
C
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F08FEY (DLATRD)
C
      END
