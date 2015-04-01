      SUBROUTINE G05GAF(SIDE,INIT,M,N,A,LDA,WK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Based on a LAPACK test routine
C
C
C     Purpose
C     =======
C
C     GO5GAF pre- or post-multiplies an M by N matrix A by a random
C     orthogonal matrix U, overwriting A. A may optionally be
C     initialized to the identity matrix before multiplying by U.
C     U is generated using the method of G.W. Stewart
C     ( SIAM J. Numer. Anal. 17, 1980, pp. 403-409 ).
C
C     Arguments
C     =========
C
C     SIDE   - CHARACTER*1
C           SIDE specifies whether A is multiplied on the left or right
C           by U.
C              SIDE = 'L'     Multiply A on the left (premultiply)
C              SIDE = 'R'     Multiply A on the right (postmultiply)
C           Unchanged on exit.
C
C     INIT   - CHARACTER*1
C           INIT specifies whether or not A should be initialized to
C           the identity matrix.
C              INIT = 'I'     Initialize A to (a section of) the
C                             identity matrix before applying U.
C              INIT = 'N'     No initialization.  Apply U to the
C                             input matrix A.
C           Unchanged on exit.
C
C     M      - INTEGER
C           Number of rows of A. Unchanged on exit.
C
C     N      - INTEGER
C           Number of columns of A. Unchanged on exit.
C
C     A      - REAL               array of DIMENSION ( LDA, N )
C           Input array. Overwritten by U*A ( if SIDE = 'L' )
C           or by A*U ( if SIDE = 'R' ) on exit.
C
C     LDA    - INTEGER
C           Leading dimension of A. Must be at least M.
C           Unchanged on exit.
C
C     WK      - REAL               array of DIMENSION (2 * MAX( M, N ) )
C           Workspace. Of length 2*M if SIDE = 'L' and of length 2*N
C           if SIDE = 'R'. Overwritten on exit.
C
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G05GAF')
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, M, N
      CHARACTER         INIT, SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,N), WK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  DUMMY, PI, S, SP
      INTEGER           I, IERROR, IROW, J, JCOL, K, LIM, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  F06EAF, F06EJF, G05CAF, G05DDF
      INTEGER           P01ABF
      EXTERNAL          F06EAF, F06EJF, G05CAF, G05DDF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06ECF, F06EDF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN, SIGN
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (M.LT.1) THEN
         WRITE (P01REC(1),FMT=99999) M
      ELSE IF (N.LT.1) THEN
         WRITE (P01REC(1),FMT=99998) N
      ELSE IF (LDA.LT.M) THEN
         WRITE (P01REC(1),FMT=99997) LDA, M
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.0) THEN
         IF (SIDE.EQ.'L' .OR. SIDE.EQ.'l') THEN
            LIM = M
         ELSE IF (SIDE.EQ.'R' .OR. SIDE.EQ.'r') THEN
            LIM = N
         ELSE
            IERROR = 2
            WRITE (P01REC(1),FMT=99996) SIDE
            GO TO 200
         END IF
C
C        Initialize A to the identity matrix if desired
C
         IF (INIT.EQ.'I' .OR. INIT.EQ.'i') THEN
            DO 40 I = 1, M
               DO 20 J = 1, N
                  A(I,J) = ZERO
   20          CONTINUE
   40       CONTINUE
            DO 60 I = 1, MIN(N,M)
               A(I,I) = ONE
   60       CONTINUE
         ELSE IF (INIT.NE.'N' .AND. INIT.NE.'n') THEN
            IERROR = 2
            WRITE (P01REC(1),FMT=99995) INIT
            GO TO 200
         END IF
         IF (LIM.EQ.1) THEN
            IERROR = 3
            WRITE (P01REC(1),FMT=99994)
            GO TO 200
         END IF
         DO 140 I = 2, LIM
            K = LIM - I + 1
            DO 80 J = K, LIM
C
C              Generate independent normal( 0, 1 ) random numbers
C
               WK(J) = G05DDF(0.0D0,1.0D0)
   80       CONTINUE
C
C           Generate a Householder transformation
C           from the random vector WK
C
            S = F06EJF(I,WK(K),1)
            S = SIGN(S,WK(K))
            WK(LIM+K) = SIGN(ONE,-S)
            PI = S*(S+WK(K))
            WK(K) = WK(K) + S
C
C           Apply Householder transformation to A
C
            IF (SIDE.EQ.'L' .OR. SIDE.EQ.'l') THEN
               DO 100 JCOL = 1, N
                  SP = -F06EAF(I,WK(K),1,A(K,JCOL),1)/PI
                  CALL F06ECF(I,SP,WK(K),1,A(K,JCOL),1)
  100          CONTINUE
            ELSE
               DO 120 IROW = 1, M
                  SP = -F06EAF(I,WK(K),1,A(IROW,K),LDA)/PI
                  CALL F06ECF(I,SP,WK(K),1,A(IROW,K),LDA)
  120          CONTINUE
            END IF
  140    CONTINUE
         WK(2*LIM) = SIGN(ONE,G05CAF(DUMMY))
C
C        Scale the matrix A by D.
C
         IF (SIDE.EQ.'L' .OR. SIDE.EQ.'l') THEN
            DO 160 IROW = 1, M
               CALL F06EDF(N,WK(LIM+IROW),A(IROW,1),LDA)
  160       CONTINUE
         ELSE
            DO 180 JCOL = 1, N
               CALL F06EDF(M,WK(LIM+JCOL),A(1,JCOL),1)
  180       CONTINUE
         END IF
      END IF
  200 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
99999 FORMAT (' ** On entry, M.lt.1: M = ',I16)
99998 FORMAT (' ** On entry, N.lt.1: N = ',I16)
99997 FORMAT (' ** On entry, LDA.lt.M: LDA = ',I16,' M = ',I16)
99996 FORMAT (' ** On entry, SIDE is not valid: SIDE = ',A1)
99995 FORMAT (' ** On entry, INIT is not valid: INIT = ',A1)
99994 FORMAT (' ** Orthogonal matrix of dimension 1 requested')
      END
