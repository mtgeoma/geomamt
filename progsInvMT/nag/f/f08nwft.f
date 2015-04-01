      SUBROUTINE F08NWF(JOB,SIDE,N,ILO,IHI,SCALE,M,V,LDV,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZGEBAK(JOB,SIDE,N,ILO,IHI,SCALE,M,V,LDV,INFO)
C
C  Purpose
C  =======
C
C  ZGEBAK forms the right or left eigenvectors of a complex general
C  matrix by backward transformation on the computed eigenvectors of the
C  matrix preprocessed by ZGEBAL.
C
C  Arguments
C  =========
C
C  JOB     (input) CHARACTER*1
C          Specifies the type of backward transformation required:
C          = 'N', do nothing, return immediately.
C          = 'P', do backward transformation for permutation.
C          = 'S', do backward transformation for balancing.
C          = 'B', do backward transformations for both permutation and
C                 balancing.
C          JOB must be the same as the JOB parameter in ZGEBAL.
C
C  SIDE    (input) CHARACTER*1
C          Specifies whether the eigenvectors given in the array V
C          are right eigenvectors or left eigenvectors.
C          = 'R', right eigenvectors
C          = 'L', left eigenvectors
C
C  N       (input) INTEGER
C          The number of rows of the matrix V.  N >= 0.
C
C  ILO     (input) INTEGER
C  IHI     (input) INTEGER
C          ILO and IHI are integers determined by ZGEBAL.
C
C  SCALE   (input) DOUBLE PRECISION array, dimension (N)
C          SCALE contains information determining the permutations
C          and/or scaling factors used by ZGEBAL.
C
C  M       (input) INTEGER
C          M is the number of columns of the matrix of eigenvectors V to
C          be back transformed.
C
C  V       (input/output) COMPLEX*16 array, dimension (LDV,M)
C          On entry, V contains the real and imaginary parts of the
C          eigenvectors to be backward transformed in its first M
C          columns.
C          On exit, V contains the real and imaginary parts of the
C          transformed eigenvectors in its first M columns.
C
C  LDV     (input) INTEGER
C          The leading dimension of the matrix V. LDV >= max(1,N).
C
C  INFO    (output) INTEGER
C          = 0:  successful exit.
C          < 0:  if INFO = -k, the k-th argument had an illegal value.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           IHI, ILO, INFO, LDV, M, N
      CHARACTER         JOB, SIDE
C     .. Array Arguments ..
      COMPLEX*16        V(LDV,*)
      DOUBLE PRECISION  SCALE(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  S
      INTEGER           I, II, K
      LOGICAL           LEFTV, RIGHTV
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, ZDSCAL, ZSWAP
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     Decode and Test the input parameters
C
      RIGHTV = (SIDE.EQ.'R' .OR. SIDE.EQ.'r')
      LEFTV = (SIDE.EQ.'L' .OR. SIDE.EQ.'l')
C
      INFO = 0
      IF ( .NOT. (JOB.EQ.'N' .OR. JOB.EQ.'n')
     *    .AND. .NOT. (JOB.EQ.'P' .OR. JOB.EQ.'p')
     *    .AND. .NOT. (JOB.EQ.'S' .OR. JOB.EQ.'s')
     *    .AND. .NOT. (JOB.EQ.'B' .OR. JOB.EQ.'b')) THEN
         INFO = -1
      ELSE IF ( .NOT. RIGHTV .AND. .NOT. LEFTV) THEN
         INFO = -2
      ELSE IF (N.LT.0) THEN
         INFO = -3
      ELSE IF (ILO.LT.1 .OR. ILO.GT.MAX(1,N)) THEN
         INFO = -4
      ELSE IF (IHI.LT.MIN(ILO,N) .OR. IHI.GT.N) THEN
         INFO = -5
      ELSE IF (M.LT.0) THEN
         INFO = -7
      ELSE IF (LDV.LT.MAX(1,N)) THEN
         INFO = -9
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08NWF/ZGEBAK',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
      IF (M.EQ.0) RETURN
      IF ((JOB.EQ.'N' .OR. JOB.EQ.'n')) RETURN
C
      IF (ILO.EQ.IHI) GO TO 60
C
C     Backward balance
C
      IF ((JOB.EQ.'S' .OR. JOB.EQ.'s') .OR. (JOB.EQ.'B' .OR. JOB.EQ.'b')
     *    ) THEN
C
         IF (RIGHTV) THEN
            DO 20 I = ILO, IHI
               S = SCALE(I)
               CALL ZDSCAL(M,S,V(I,1),LDV)
   20       CONTINUE
         END IF
C
         IF (LEFTV) THEN
            DO 40 I = ILO, IHI
               S = ONE/SCALE(I)
               CALL ZDSCAL(M,S,V(I,1),LDV)
   40       CONTINUE
         END IF
C
      END IF
C
C     Backward permutation
C
C     For  I = ILO-1 step -1 until 1,
C              IHI+1 step 1 until N do --
C
   60 CONTINUE
      IF ((JOB.EQ.'P' .OR. JOB.EQ.'p') .OR. (JOB.EQ.'B' .OR. JOB.EQ.'b')
     *    ) THEN
         IF (RIGHTV) THEN
            DO 80 II = 1, N
               I = II
               IF (I.GE.ILO .AND. I.LE.IHI) GO TO 80
               IF (I.LT.ILO) I = ILO - II
               K = SCALE(I)
               IF (K.EQ.I) GO TO 80
               CALL ZSWAP(M,V(I,1),LDV,V(K,1),LDV)
   80       CONTINUE
         END IF
C
         IF (LEFTV) THEN
            DO 100 II = 1, N
               I = II
               IF (I.GE.ILO .AND. I.LE.IHI) GO TO 100
               IF (I.LT.ILO) I = ILO - II
               K = SCALE(I)
               IF (K.EQ.I) GO TO 100
               CALL ZSWAP(M,V(I,1),LDV,V(K,1),LDV)
  100       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F08NWF (ZGEBAK)
C
      END
