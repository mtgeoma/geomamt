      SUBROUTINE Y90CDF(SIDE,INIT,M,N,A,LDA,X,D,SEED)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C  -- LAPACK test routine --
C     Argonne National Lab and Courant Institute
C     November 14, 1988
C
C
C  Purpose
C  =======
C
C     Y90CDF pre- or post-multiplies an M by N matrix A by a random
C     unitary matrix U, overwriting A. A may optionally be
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
C  A      - COMPLEX            array of DIMENSION ( LDA, N )
C           Input array. Overwritten by U*A ( if SIDE = 'L' )
C           or by A*U ( if SIDE = 'R' ) on exit.
C
C  LDA    - INTEGER
C           Leading dimension of A. Must be at least MAX ( 1, M ).
C           Unchanged on exit.
C
C  X      - COMPLEX            array of DIMENSION ( MAX( M, N ) )
C           Workspace. Of length M if SIDE = 'L' and of length N
C           if SIDE = 'R'. Overwritten on exit.
C
C  D      - COMPLEX            array of DIMENSION ( MAX( M, N ) )
C           Workspace. Of length M if SIDE = 'L' and of length N
C           if SIDE = 'R'. Overwritten on exit.
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      COMPLEX*16        CZERO, CONE
      PARAMETER         (CZERO=0.0D0,CONE=1.0D0)
C     .. Scalar Arguments ..
      INTEGER           LDA, M, N
      CHARACTER*1       INIT, SIDE
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), D(*), X(*)
      INTEGER           SEED(4)
C     .. Local Scalars ..
      COMPLEX*16        CSIGN, SP, XNORMS
      DOUBLE PRECISION  FACTOR, XABS, XNORM
      INTEGER           I, IROW, IXFRM, J, JCOL, KBEG, NXFRM
C     .. External Functions ..
      COMPLEX*16        F06GBF, Y90EBF
      DOUBLE PRECISION  F06JJF
      LOGICAL           Y90WAF
      EXTERNAL          F06GBF, Y90EBF, F06JJF, Y90WAF
C     .. External Subroutines ..
      EXTERNAL          F06GCF, F06GDF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DCONJG, MIN
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C
      IF (N.LT.1 .OR. M.LT.1) RETURN
      IF (Y90WAF(SIDE,'L')) THEN
         NXFRM = M
      ELSE
         NXFRM = N
      END IF
C
C     Initialize A to the identity matrix if desired
C
      IF (Y90WAF(INIT,'I')) THEN
C
         DO 40 I = 1, M
            DO 20 J = 1, N
               A(I,J) = ZERO
   20       CONTINUE
   40    CONTINUE
C
         DO 60 I = 1, MIN(N,M)
            A(I,I) = CONE
   60    CONTINUE
      END IF
C
C     Return if no rotation possible
C
      IF (NXFRM.EQ.1) RETURN
C
C-----------------------------------------------------------------------
C
C
C        2)      Compute Rotation by computing Householder Transformatio
C                H(2), H(3), ..., H(n).  Note that the order in which
C                they are computed is irrelevant.
C
C
      DO 80 J = 1, NXFRM
         X(J) = ZERO
   80 CONTINUE
C
C
      DO 160 IXFRM = 2, NXFRM
         KBEG = NXFRM - IXFRM + 1
C
C           Generate independent normal( 0, 1 ) random numbers
C
         DO 100 J = KBEG, NXFRM
            X(J) = Y90EBF(4,SEED)
  100    CONTINUE
C
C        Generate a Householder transformation from the random vector X
C
         XNORM = F06JJF(IXFRM,X(KBEG),1)
         XABS = ABS(X(KBEG))
         IF (XABS.NE.CZERO) THEN
            CSIGN = X(KBEG)/XABS
         ELSE
            CSIGN = CONE
         END IF
         XNORMS = CSIGN*XNORM
         D(KBEG) = -CSIGN
         FACTOR = XNORM*(XNORM+XABS)
         X(KBEG) = X(KBEG) + XNORMS
C
C        Apply Householder transformation to A
C
         IF (Y90WAF(SIDE,'L')) THEN
C
            DO 120 JCOL = 1, N
               SP = -F06GBF(IXFRM,X(KBEG),1,A(KBEG,JCOL),1)/FACTOR
               CALL F06GCF(IXFRM,SP,X(KBEG),1,A(KBEG,JCOL),1)
  120       CONTINUE
C
         ELSE
C
            DO 140 IROW = 1, M
               SP = -F06GBF(IXFRM,X(KBEG),1,A(IROW,KBEG),LDA)/FACTOR
               CALL F06GCF(IXFRM,SP,X(KBEG),1,A(IROW,KBEG),LDA)
  140       CONTINUE
C
         END IF
  160 CONTINUE
C
      X(1) = Y90EBF(4,SEED)
      XABS = ABS(X(1))
      IF (XABS.NE.CZERO) THEN
         CSIGN = X(KBEG)/XABS
      ELSE
         CSIGN = CONE
      END IF
      D(NXFRM) = CSIGN
C
C-----------------------------------------------------------------------
C
C     Scale the matrix A by D.
C
      IF (Y90WAF(SIDE,'L')) THEN
         DO 180 IROW = 1, M
            CALL F06GDF(N,DCONJG(D(IROW)),A(IROW,1),LDA)
  180    CONTINUE
      ELSE
         DO 200 JCOL = 1, N
            CALL F06GDF(M,D(JCOL),A(1,JCOL),1)
  200    CONTINUE
      END IF
      RETURN
C
C     End of Y90CDF
      END
