      SUBROUTINE F08NTF(N,ILO,IHI,A,LDA,TAU,WORK,LWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZUNGHR(N,ILO,IHI,A,LDA,TAU,WORK,LWORK,INFO)
C
C  Purpose
C  =======
C
C  ZUNGHR generates a complex unitary matrix Q which is defined as the
C  product of ihi-ilo elementary reflectors of order n, as returned by
C  ZGEHRD:
C
C  Q = H(ilo) H(ilo+1) . . . H(ihi-1).
C
C  Arguments
C  =========
C
C  N       (input) INTEGER
C          The order of the matrix Q. N >= 0.
C
C  ILO     (input) INTEGER
C  IHI     (input) INTEGER
C          ILO and IHI must have the same values as in the previous call
C          of ZGEHRD. Q is equal to the unit matrix except in the
C          submatrix Q(ilo+1:ihi,ilo+1:ihi).  If N > 0,
C          1 <= ILO <= IHI <= N; otherwise ILO = 1 and IHI = N.
C
C  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C          On entry, the vectors which define the elementary reflectors,
C          as returned by ZGEHRD.
C          On exit, the n by n unitary matrix Q.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A. LDA >= max(1,N).
C
C  TAU     (input) COMPLEX*16 array, dimension (N-1)
C          TAU(i) must contain the scalar factor of the elementary
C          reflector H(i), as returned by ZGEHRD.
C
C  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
C          On exit, if INFO = 0, WORK(1) returns the minimum value of
C          LWORK required to use the optimal blocksize.
C
C  LWORK   (input) INTEGER
C          The dimension of the array WORK. LWORK >= IHI-ILO.
C          For optimum performance LWORK should be at least
C          (IHI-ILO)*NB, where NB is the optimal blocksize.
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument had an illegal value
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           IHI, ILO, INFO, LDA, LWORK, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), TAU(*), WORK(LWORK)
C     .. Local Scalars ..
      INTEGER           I, IINFO, J, NH
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, ZUNGQR
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     Test the input arguments
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
      ELSE IF (LWORK.LT.MAX(1,IHI-ILO)) THEN
         INFO = -8
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08NTF/ZUNGHR',-INFO)
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
C     Shift the vectors which define the elementary reflectors one
C     column to the right, and set the first ilo and the last n-ihi
C     rows and columns to those of the unit matrix
C
      DO 80 J = IHI, ILO + 1, -1
         DO 20 I = 1, J - 1
            A(I,J) = ZERO
   20    CONTINUE
         DO 40 I = J + 1, IHI
            A(I,J) = A(I,J-1)
   40    CONTINUE
         DO 60 I = IHI + 1, N
            A(I,J) = ZERO
   60    CONTINUE
   80 CONTINUE
      DO 120 J = 1, ILO
         DO 100 I = 1, N
            A(I,J) = ZERO
  100    CONTINUE
         A(J,J) = ONE
  120 CONTINUE
      DO 160 J = IHI + 1, N
         DO 140 I = 1, N
            A(I,J) = ZERO
  140    CONTINUE
         A(J,J) = ONE
  160 CONTINUE
C
      NH = IHI - ILO
      IF (NH.GT.0) THEN
C
C        Generate Q(ilo+1:ihi,ilo+1:ihi)
C
         CALL ZUNGQR(NH,NH,NH,A(ILO+1,ILO+1),LDA,TAU(ILO),WORK,LWORK,
     *               IINFO)
      ELSE
         WORK(1) = 1
      END IF
      RETURN
C
C     End of F08NTF (ZUNGHR)
C
      END
