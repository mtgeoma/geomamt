      SUBROUTINE Y90RTX(N,NJORD,JORD,EIGR,EIGI,A,LDA,X,LDX,Y,LDY,WK1,
     *                  LDWK1)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C-----------------------------------------------------------------------
C
C         ====================================
C         *  Y90RTX :  Auxiliary for Y90RTF  *
C         ====================================
C
C
C     DESCRIPTION
C     ===========
C
C     Y90RTX ...
C
C
C     DEPENDENCIES
C     ============
C
C     None
C
C
C     ARGUMENTS
C     =========
C
C     N      :  Integer, input.
C               Order of the matrices to be generated.
C
C     NJORD  :  Integer, input.
C               Number of Jordan blocks.
C
C     JORD   :  Integer array of DIMENSION (*), input.
C               The absolute value of JORD(i) denotes the size of the
C               i-th Jordan block.  For complex eigenpairs, JORD(i)
C               is negative and only one of the pair needs to be
C               specified, and only one Jordan block.
C
C     EIGR   :  Real array of DIMENSION (*), input/output.
C               EIGR(i) contains the real part of the i-th eigenvalue.
C
C     EIGI   :  Real array of DIMENSION (*), input/output.
C               EIGR(i) contains the imaginary part of the i-th
C               eigenvalue.
C
C     A      :  Real array of DIMENSION (LDA,*), output.
C               The matrix to be generated.  It must have at least N
C               columns.
C
C     LDA    :  Integer, input.
C               Leading dimension of the array A.  It must be LDA >= N.
C
C     X      :  Real array of DIMENSION (LDX,*), input.
C               Contains the matrix X. It must have at least N columns.
C
C     LDX    :  Integer, input.
C               Leading dimension of the array X.  It must be LDX >= N.
C
C     Y      :  Real array of DIMENSION (LDY,*), input.
C               Contains the transpose of the matrix X**(-1). It must
C               have at least N columns.
C
C     LDY    :  Integer, input.
C               Leading dimension of the array Y. It must be LDY >= N.
C
C     WK1    :  Real array of DIMENSION (LDWK1,*), workspace.
C               It must have at least N columns.
C
C     LDWK1  :  Integer, input.
C               Leading dimension of WK1.  It must be LDWK1 >= N.
C
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      INTEGER           LDA, LDWK1, LDX, LDY, N, NJORD
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), EIGI(*), EIGR(*), WK1(LDWK1,*),
     *                  X(LDX,*), Y(LDY,*)
      INTEGER           JORD(*)
C     .. Local Scalars ..
      INTEGER           I, J, K
      LOGICAL           PAIR
C     .. External Subroutines ..
      EXTERNAL          F06ECF, F06EDF, F06YAF, Y90SHF
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Carry out the computation
C
C-----------------------------------------------------------------------
      CALL Y90SHF('N',N,N,X,LDX,WK1,LDWK1)
      PAIR = .FALSE.
      DO 20 I = 1, N
         IF (PAIR) THEN
            PAIR = .FALSE.
         ELSE
            CALL F06EDF(N,EIGR(I),WK1(1,I),1)
            IF (EIGI(I).NE.ZERO) THEN
               PAIR = .TRUE.
               CALL F06EDF(N,EIGR(I),WK1(1,I+1),1)
               CALL F06ECF(N,-EIGI(I),X(1,I+1),1,WK1(1,I),1)
               CALL F06ECF(N,EIGI(I),X(1,I),1,WK1(1,I+1),1)
            END IF
         END IF
   20 CONTINUE
C
      K = 1
      DO 80 I = 1, NJORD
         IF (JORD(I).GT.0) THEN
            DO 40 J = K, K + JORD(I) - 2
               CALL F06ECF(N,ONE,X(1,J),1,WK1(1,J+1),1)
   40       CONTINUE
            K = K + JORD(I)
         ELSE
            DO 60 J = K, K + 2*(-JORD(I)-1) - 1
               CALL F06ECF(N,ONE,X(1,J),1,WK1(1,J+2),1)
   60       CONTINUE
            K = K - 2*JORD(I)
         END IF
   80 CONTINUE
C
      CALL F06YAF('N','T',N,N,N,ONE,WK1,LDWK1,Y,LDY,ZERO,A,LDA)
C-----------------------------------------------------------------------
C
C     End of subroutine Y90RTX
C
C-----------------------------------------------------------------------
      RETURN
      END
