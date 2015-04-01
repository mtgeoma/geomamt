      SUBROUTINE F08NSF(N,ILO,IHI,A,LDA,TAU,WORK,LWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZGEHRD(N,ILO,IHI,A,LDA,TAU,WORK,LWORK,INFO)
C
C  Purpose
C  =======
C
C  ZGEHRD reduces a complex general matrix A to upper Hessenberg form H
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
C          Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to
C          zero.
C
C  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
C          On exit, if INFO = 0, WORK(1) returns the minimum value of
C          LWORK required to use the optimum blocksize.
C
C  LWORK   (input) INTEGER
C          The length of the array WORK.  LWORK >= max(1,N).
C          For optimum performance LWORK should be at least N*NB,
C          where NB is the optimal blocksize.
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
      INTEGER           NBMAX, LDT
      PARAMETER         (NBMAX=64,LDT=NBMAX+1)
      COMPLEX*16        ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           IHI, ILO, INFO, LDA, LWORK, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), TAU(*), WORK(LWORK)
C     .. Local Scalars ..
      COMPLEX*16        EI
      INTEGER           I, IB, IINFO, IWS, LDWORK, NB, NBMIN, NH, NX
C     .. Local Arrays ..
      COMPLEX*16        T(LDT,NBMAX)
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07ZAZ, F08ASY, F08NSY, F08NSZ, ZGEMM
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
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
      ELSE IF (LWORK.LT.MAX(1,N)) THEN
         INFO = -8
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08NSF/ZGEHRD',-INFO)
         RETURN
      END IF
C
C     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero
C
      DO 20 I = 1, ILO - 1
         TAU(I) = ZERO
   20 CONTINUE
      DO 40 I = MAX(1,IHI), N - 1
         TAU(I) = ZERO
   40 CONTINUE
C
C     Quick return if possible
C
      NH = IHI - ILO + 1
      IF (NH.LE.1) THEN
         WORK(1) = 1
         RETURN
      END IF
C
C     Determine the block size.
C
      CALL F07ZAZ(1,'F08NSF',NB,0)
      NB = MIN(NBMAX,NB)
      NBMIN = 2
      IWS = 1
      IF (NB.GT.1 .AND. NB.LT.NH) THEN
C
C        Determine when to cross over from blocked to unblocked code
C        (last block is always handled by unblocked code).
C
         CALL F07ZAZ(3,'F08NSF',NX,0)
         NX = MAX(NB,NX)
         IF (NX.LT.NH) THEN
C
C           Determine if workspace is large enough for blocked code.
C
            IWS = N*NB
            IF (LWORK.LT.IWS) THEN
C
C              Not enough workspace to use optimal NB:  determine the
C              minimum value of NB, and reduce NB or force use of
C              unblocked code.
C
               CALL F07ZAZ(2,'F08NSF',NBMIN,0)
               NBMIN = MAX(2,NBMIN)
               IF (LWORK.GE.N*NBMIN) THEN
                  NB = LWORK/N
               ELSE
                  NB = 1
               END IF
            END IF
         END IF
      END IF
      LDWORK = N
C
      IF (NB.LT.NBMIN .OR. NB.GE.NH) THEN
C
C        Use unblocked code below
C
         I = ILO
C
      ELSE
C
C        Use blocked code
C
         DO 60 I = ILO, IHI - 1 - NX, NB
            IB = MIN(NB,IHI-I)
C
C           Reduce columns i:i+ib-1 to Hessenberg form, returning the
C           matrices V and T of the block reflector H = I - V*T*V'
C           which performs the reduction, and also the matrix Y = A*V*T
C
            CALL F08NSY(IHI,I,IB,A(1,I),LDA,TAU(I),T,LDT,WORK,LDWORK)
C
C           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
C           right, computing  A := A - Y * V'. V(i+ib,ib-1) must be set
C           to 1.
C
            EI = A(I+IB,I+IB-1)
            A(I+IB,I+IB-1) = ONE
            CALL ZGEMM('No transpose','Conjugate transpose',IHI,
     *                 IHI-I-IB+1,IB,-ONE,WORK,LDWORK,A(I+IB,I),LDA,ONE,
     *                 A(1,I+IB),LDA)
            A(I+IB,I+IB-1) = EI
C
C           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
C           left
C
            CALL F08ASY('Left','Conjugate transpose','Forward',
     *                  'Columnwise',IHI-I,N-I-IB+1,IB,A(I+1,I),LDA,T,
     *                  LDT,A(I+1,I+IB),LDA,WORK,LDWORK)
   60    CONTINUE
      END IF
C
C     Use unblocked code to reduce the rest of the matrix
C
      CALL F08NSZ(N,I,IHI,A,LDA,TAU,WORK,IINFO)
      WORK(1) = IWS
C
      RETURN
C
C     End of F08NSF (ZGEHRD)
C
      END
