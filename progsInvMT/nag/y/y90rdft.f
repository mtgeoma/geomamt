      SUBROUTINE Y90RDF(SIDE,INIT,M,N,A,LDA,X,D,SEED)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C  -- LAPACK test routine --
C     Argonne National Lab and Courant Institute
C     October 11, 1988
C
C
C  Purpose
C  =======
C
C     Y90RDF pre- or post-multiplies an M by N matrix A by a random
C     orthogonal matrix U, overwriting A. A may optionally be
C     initialized to the identity matrix before multiplying by U.
C     U is generated using the method of G.W. Stewart
C     ( SIAM J. Numer. Anal. 17, 1980, pp. 403-409 ).
C
C  Arguments
C  =========
C
C  SIDE   - CHARACTER*1
C           SIDE specifies whether A is multiplied on the left or right
C           by U.
C              SIDE = 'L'     Multiply A on the left (premultiply)
C              SIDE = 'U'     Multiply A on the right (postmultiply)
C           Unchanged on exit.
C
C  INIT   - CHARACTER*1
C           INIT specifies whether or not A should be initialized to
C           the identity matrix.
C              INIT = 'I'     Initialize A to (a section of) the
C                             identity matrix before applying U.
C              INIT = 'N'     No initialization.  Apply U to the
C                             input matrix A.
C           Unchanged on exit.
C
C  M      - INTEGER
C           Number of rows of A. Unchanged on exit.
C
C  N      - INTEGER
C           Number of columns of A. Unchanged on exit.
C
C  A      - REAL               array of DIMENSION ( LDA, N )
C           Input array. Overwritten by U*A ( if SIDE = 'L' )
C           or by A*U ( if SIDE = 'R' ) on exit.
C
C  LDA    - INTEGER
C           Leading dimension of A. Must be at least MAX ( 1, M ).
C           Unchanged on exit.
C
C  X      - REAL               array of DIMENSION ( MAX( M, N ) )
C           Workspace. Of length M if SIDE = 'L' and of length N
C           if SIDE = 'R'. Overwritten on exit.
C
C  D      - REAL               array of DIMENSION ( MAX( M, N ) )
C           Workspace. Of length M if SIDE = 'L' and of length N
C           if SIDE = 'R'. Overwritten on exit.
C
C
C
C
C
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           LDA, M, N
      CHARACTER*1       INIT, SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), D(*), X(*)
      INTEGER           SEED(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  PI, S, SP
      INTEGER           I, IROW, J, JCOL, K, LIM
C     .. External Functions ..
      DOUBLE PRECISION  F06EAF, F06EJF, Y90TBF
      LOGICAL           Y90WAF
      EXTERNAL          F06EAF, F06EJF, Y90TBF, Y90WAF
C     .. External Subroutines ..
      EXTERNAL          F06ECF, F06EDF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN, SIGN
C     .. Executable Statements ..
C
      IF (N.LT.1 .OR. M.LT.1) RETURN
      IF (Y90WAF(SIDE,'L')) THEN
         LIM = M
      ELSE
         LIM = N
      END IF
C
C     Initialize A to the identity matrix if desired
C
      IF (Y90WAF(INIT,'I')) THEN
         DO 40 I = 1, M
            DO 20 J = 1, N
               A(I,J) = ZERO
   20       CONTINUE
   40    CONTINUE
         DO 60 I = 1, MIN(N,M)
            A(I,I) = ONE
   60    CONTINUE
      END IF
      IF (LIM.EQ.1) RETURN
      DO 140 I = 2, LIM
         K = LIM - I + 1
         DO 80 J = K, LIM
C
C           Generate independent normal( 0, 1 ) random numbers
C
            X(J) = Y90TBF(3,SEED)
   80    CONTINUE
C
C        Generate a Householder transformation from the random vector X
C
         S = F06EJF(I,X(K),1)
         S = SIGN(S,X(K))
         D(K) = SIGN(ONE,-S)
         PI = S*(S+X(K))
         X(K) = X(K) + S
C
C        Apply Householder transformation to A
C
         IF (Y90WAF(SIDE,'L')) THEN
            DO 100 JCOL = 1, N
               SP = -F06EAF(I,X(K),1,A(K,JCOL),1)/PI
               CALL F06ECF(I,SP,X(K),1,A(K,JCOL),1)
  100       CONTINUE
         ELSE
            DO 120 IROW = 1, M
               SP = -F06EAF(I,X(K),1,A(IROW,K),LDA)/PI
               CALL F06ECF(I,SP,X(K),1,A(IROW,K),LDA)
  120       CONTINUE
         END IF
  140 CONTINUE
      D(LIM) = SIGN(ONE,Y90TBF(3,SEED))
C
C     Scale the matrix A by D.
C
      IF (Y90WAF(SIDE,'L')) THEN
         DO 160 IROW = 1, M
            CALL F06EDF(N,D(IROW),A(IROW,1),LDA)
  160    CONTINUE
      ELSE
         DO 180 JCOL = 1, N
            CALL F06EDF(M,D(JCOL),A(1,JCOL),1)
  180    CONTINUE
      END IF
      RETURN
C
C     End of Y90RDF
      END
