      SUBROUTINE F08JEV(TYPE,KL,KU,CFROM,CTO,M,N,A,LDA,INFO)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C  Purpose
C  =======
C
C  DLASCL multiplies the M by N real matrix A by the real scalar
C  CTO/CFROM.  This is done without over/underflow as long as the final
C  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
C  A may be full, upper triangular, lower triangular, upper Hessenberg,
C  or banded.
C
C  Arguments
C  =========
C
C  TYPE    (input) CHARACTER*1
C          TYPE indices the storage type of the input matrix.
C          = 'G':  A is a full matrix.
C          = 'L':  A is a lower triangular matrix.
C          = 'U':  A is an upper triangular matrix.
C          = 'H':  A is an upper Hessenberg matrix.
C          = 'B':  A is a symmetric band matrix with lower bandwidth KL
C                  and upper bandwidth KU and with the only the lower
C                  half stored.
C          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
C                  and upper bandwidth KU and with the only the upper
C                  half stored.
C          = 'Z':  A is a band matrix with lower bandwidth KL and upper
C                  bandwidth KU.
C
C  KL      (input) INTEGER
C          The lower bandwidth of A.  Referenced only if TYPE = 'B',
C          'Q' or 'Z'.
C
C  KU      (input) INTEGER
C          The upper bandwidth of A.  Referenced only if TYPE = 'B',
C          'Q' or 'Z'.
C
C  CFROM   (input) DOUBLE PRECISION
C  CTO     (input) DOUBLE PRECISION
C          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
C          without over/underflow if the final result CTO*A(I,J)/CFROM
C          can be represented without over/underflow.  CFROM must be
C          nonzero.
C
C  M       (input) INTEGER
C          The number of rows of the matrix A.  M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix A.  N >= 0.
C
C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,M)
C          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
C          storage type.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,M).
C
C  INFO    (output) INTEGER
C          0  - successful exit
C          <0 - if INFO = -i, the i-th argument had an illegal value.
C
C -- LAPACK auxiliary routine (version 2.0) (adapted for NAG library) --
C    Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C    Courant Institute, Argonne National Lab, and Rice University
C    February 29, 1992
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CFROM, CTO
      INTEGER           INFO, KL, KU, LDA, M, N
      CHARACTER         TYPE
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
      INTEGER           I, ITYPE, J, K1, K2, K3, K4
      LOGICAL           DONE
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      EXTERNAL          X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
C
      IF (TYPE.EQ.'G' .OR. TYPE.EQ.'g') THEN
         ITYPE = 0
      ELSE IF (TYPE.EQ.'L' .OR. TYPE.EQ.'l') THEN
         ITYPE = 1
      ELSE IF (TYPE.EQ.'U' .OR. TYPE.EQ.'u') THEN
         ITYPE = 2
      ELSE IF (TYPE.EQ.'H' .OR. TYPE.EQ.'h') THEN
         ITYPE = 3
      ELSE IF (TYPE.EQ.'B' .OR. TYPE.EQ.'b') THEN
         ITYPE = 4
      ELSE IF (TYPE.EQ.'Q' .OR. TYPE.EQ.'q') THEN
         ITYPE = 5
      ELSE IF (TYPE.EQ.'Z' .OR. TYPE.EQ.'z') THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF
C
      IF (ITYPE.EQ.-1) THEN
         INFO = -1
      ELSE IF (CFROM.EQ.ZERO) THEN
         INFO = -4
      ELSE IF (M.LT.0) THEN
         INFO = -6
      ELSE IF (N.LT.0 .OR. (ITYPE.EQ.4 .AND. N.NE.M)
     *         .OR. (ITYPE.EQ.5 .AND. N.NE.M)) THEN
         INFO = -7
      ELSE IF (ITYPE.LE.3 .AND. LDA.LT.MAX(1,M)) THEN
         INFO = -9
      ELSE IF (ITYPE.GE.4) THEN
         IF (KL.LT.0 .OR. KL.GT.MAX(M-1,0)) THEN
            INFO = -2
         ELSE IF (KU.LT.0 .OR. KU.GT.MAX(N-1,0)
     *            .OR. ((ITYPE.EQ.4 .OR. ITYPE.EQ.5) .AND. KL.NE.KU))
     *            THEN
            INFO = -3
         ELSE IF ((ITYPE.EQ.4 .AND. LDA.LT.KL+1)
     *            .OR. (ITYPE.EQ.5 .AND. LDA.LT.KU+1)
     *            .OR. (ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1)) THEN
            INFO = -9
         END IF
      END IF
C
      IF (INFO.NE.0) THEN
C         CALL XERBLA('DLASCL',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0 .OR. M.EQ.0) RETURN
C
C     Get machine parameters
C
      SMLNUM = X02AMF()
      BIGNUM = ONE/SMLNUM
C
      CFROMC = CFROM
      CTOC = CTO
C
   20 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      CTO1 = CTOC/BIGNUM
      IF (ABS(CFROM1).GT.ABS(CTOC) .AND. CTOC.NE.ZERO) THEN
         MUL = SMLNUM
         DONE = .FALSE.
         CFROMC = CFROM1
      ELSE IF (ABS(CTO1).GT.ABS(CFROMC)) THEN
         MUL = BIGNUM
         DONE = .FALSE.
         CTOC = CTO1
      ELSE
         MUL = CTOC/CFROMC
         DONE = .TRUE.
      END IF
C
      IF (ITYPE.EQ.0) THEN
C
C        Full matrix
C
         DO 60 J = 1, N
            DO 40 I = 1, M
               A(I,J) = A(I,J)*MUL
   40       CONTINUE
   60    CONTINUE
C
      ELSE IF (ITYPE.EQ.1) THEN
C
C        Lower triangular matrix
C
         DO 100 J = 1, N
            DO 80 I = J, M
               A(I,J) = A(I,J)*MUL
   80       CONTINUE
  100    CONTINUE
C
      ELSE IF (ITYPE.EQ.2) THEN
C
C        Upper triangular matrix
C
         DO 140 J = 1, N
            DO 120 I = 1, MIN(J,M)
               A(I,J) = A(I,J)*MUL
  120       CONTINUE
  140    CONTINUE
C
      ELSE IF (ITYPE.EQ.3) THEN
C
C        Upper Hessenberg matrix
C
         DO 180 J = 1, N
            DO 160 I = 1, MIN(J+1,M)
               A(I,J) = A(I,J)*MUL
  160       CONTINUE
  180    CONTINUE
C
      ELSE IF (ITYPE.EQ.4) THEN
C
C        Lower half of a symmetric band matrix
C
         K3 = KL + 1
         K4 = N + 1
         DO 220 J = 1, N
            DO 200 I = 1, MIN(K3,K4-J)
               A(I,J) = A(I,J)*MUL
  200       CONTINUE
  220    CONTINUE
C
      ELSE IF (ITYPE.EQ.5) THEN
C
C        Upper half of a symmetric band matrix
C
         K1 = KU + 2
         K3 = KU + 1
         DO 260 J = 1, N
            DO 240 I = MAX(K1-J,1), K3
               A(I,J) = A(I,J)*MUL
  240       CONTINUE
  260    CONTINUE
C
      ELSE IF (ITYPE.EQ.6) THEN
C
C        Band matrix
C
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 300 J = 1, N
            DO 280 I = MAX(K1-J,K2), MIN(K3,K4-J)
               A(I,J) = A(I,J)*MUL
  280       CONTINUE
  300    CONTINUE
C
      END IF
C
      IF ( .NOT. DONE) GO TO 20
C
      RETURN
C
C     End of F08JEV (DLASCL)
C
      END
