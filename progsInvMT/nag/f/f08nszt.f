      SUBROUTINE F08NSZ(N,ILO,IHI,A,LDA,TAU,WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZGEHD2(N,ILO,IHI,A,LDA,TAU,WORK,INFO)
C
C  Purpose
C  =======
C
C  ZGEHD2 reduces a complex general matrix A to upper Hessenberg form H
C  by a unitary similarity transformation:  Q' * A * Q = H .
C
C  Arguments
C  =========
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  ILO     (input) INTEGER
C  IHI     (input) INTEGER
C          It is assumed that A is already upper triangular in rows
C          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
C          set by a previous call to ZGEBAL; otherwise they should be
C          set to 1 and N respectively. See Further Details. If N > 0,
C          1 <= ILO <= IHI <= N; if N = 0, ILO = 1 and IHI = 0.
C
C  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C          On entry, the n by n general matrix to be reduced.
C          On exit, the upper triangle and the first subdiagonal of A
C          are overwritten with the upper Hessenberg matrix H, and the
C          elements below the first subdiagonal, with the array TAU,
C          represent the unitary matrix Q as a product of elementary
C          reflectors. See Further Details.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  TAU     (output) COMPLEX*16 array, dimension (N-1)
C          The scalar factors of the elementary reflectors (see Further
C          Details).
C
C  WORK    (workspace) COMPLEX*16 array, dimension (N)
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0:  if INFO = -i, the i-th argument had an illegal value.
C
C  Further Details
C  ===============
C
C  The matrix Q is represented as a product of (ihi-ilo) elementary
C  reflectors
C
C     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
C
C  Each H(i) has the form
C
C     H(i) = I - tau * v * v'
C
C  where tau is a complex scalar, and v is a complex vector with
C  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
C  exit in A(i+2:ihi,i), and tau in TAU(i).
C
C  The contents of A are illustrated by the following example, with
C  n = 7, ilo = 2 and ihi = 6:
C
C  on entry                         on exit
C
C  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
C  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
C  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
C  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
C  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
C  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
C  (                         a )    (                          a )
C
C  where a denotes an element of the original matrix A, h denotes a
C  modified element of the upper Hessenberg matrix H, and vi denotes an
C  element of the vector defining H(i).
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           IHI, ILO, INFO, LDA, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), TAU(*), WORK(*)
C     .. Local Scalars ..
      COMPLEX*16        ALPHA
      INTEGER           I
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F08ASV, F08ASW
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG, MAX, MIN
C     .. Executable Statements ..
C
C     Test the input parameters
C
      INFO = 0
      IF (N.LT.0) THEN
         INFO = -1
      ELSE IF (ILO.LT.1 .OR. ILO.GT.MAX(1,N)) THEN
         INFO = -2
      ELSE IF (IHI.LT.MIN(ILO,N) .OR. IHI.GT.N) THEN
         INFO = -3
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -5
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08NSZ/ZGEHD2',-INFO)
         RETURN
      END IF
C
      DO 20 I = ILO, IHI - 1
C
C        Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
C
         ALPHA = A(I+1,I)
         CALL F08ASV(IHI-I,ALPHA,A(MIN(I+2,N),I),1,TAU(I))
         A(I+1,I) = ONE
C
C        Apply H(i) to A(1:ihi,i+1:ihi) from the right
C
         CALL F08ASW('Right',IHI,IHI-I,A(I+1,I),1,TAU(I),A(1,I+1),LDA,
     *               WORK)
C
C        Apply H(i)' to A(i+1:ihi,i+1:n) from the left
C
         CALL F08ASW('Left',IHI-I,N-I,A(I+1,I),1,DCONJG(TAU(I)),
     *               A(I+1,I+1),LDA,WORK)
C
         A(I+1,I) = ALPHA
   20 CONTINUE
C
      RETURN
C
C     End of F08NSZ (ZGEHD2)
C
      END
