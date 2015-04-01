      SUBROUTINE Y90RQF(JOB,N,D,T0,T,Q,LDQ,ERR,WORK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         =========================================
C         *  Y90RQF :  Toeplitz Matrix Generator  *
C         =========================================
C
C
C     -- Written on 10-October-1990.
C     Sven Hammarling, Nag Ltd.
C
C
C     1. Purpose
C     ==========
C
C     Y90RQF generates an  n by n symmetric Toeplitz matrix  T  with
C     given eigenvalues. The eigenvectors of T are also returned.
C
C     2. Description
C     ==============
C
C     Y90RQF  uses  the  iterative method  described  in [1]  to
C     generate a symmetric Toeplitz matrix
C
C     T = ( t(0)    t(1)    t(2)    ...  t(n-1) ),
C         ( t(1)    t(0)    t(1)    ...  t(n-2) )
C         ( t(2)    t(1)    t(0)    ...  t(n-3) )
C         (   :       :       :            :    )
C         ( t(n-1)  t(n-2)  t(n-3)  ...  t(0)   )
C
C     with given eigenvalues.
C
C     3. References
C     =============
C
C     [1] Laurie, D. P.
C      A  numerical  approach  to  the  inverse  Toeplitz  eigenproblem.
C      SIAM J. Sci. Stat Comput., 9, 401-405, 1988.
C
C     4. Parameters
C     =============
C
C  1: JOB    - INTEGER.                                            Input
C
C              On entry: defines the task to be carried out by Y90RQF.
C              JOB = 1  :  initialize the Toeplitz matrix as in [1],
C                          using a larger tolerance and up to 5*n
C                          iterations.
C              JOB = 2  :  use a previously calculated Toeplitz matrix
C                          as a starting point,using a larger tolerance
C                          and up to 5*n iterations.
C              JOB = 3  :  use a previously calculated Toeplitz matrix
C                          as a starting point.
C              JOB = 4  :  initialize the Toeplitz matrix as in [1].
C
C  2: N      - INTEGER.                                            Input
C
C              On entry: n, the order of the required Toeplitz matrix.
C
C              Constraint: N .ge. 1.
C
C  3: D( N ) - REAL array.                                  Input/Output
C
C              On entry: the required eigenvalues of the Toeplitz matrix
C              T.
C
C              On exit:  the required eigenvalues arranged  in ascending
C              order.
C
C  4: T0     - REAL.                                              Output
C
C              On exit: the  value  of  the  diagonal  elements  of  the
C              required Toeplitz matrix.
C
C  5: T( * ) - REAL array.                                        Output
C              Note: the length of  T must be at least  max( 1, N - 1 ).
C
C              On exit: the (n-1) parameters  t(1), t(2), ..., t(n-1) of
C              the required Toeplitz matrix.
C
C  6: Q( LDQ, N ) - REAL array.                                   Output
C
C               On exit: the eigenvector matrix of the required Toeplitz
C               matrix T.  Q is orthogonal and such that
C
C                  T = Q*diag(D)*Q'.
C
C  7: LDQ    - INTEGER.                                            Input
C
C              On entry: the first dimension of the array  Q as declared
C              in the (sub)program from which TOEPLZ is called.
C
C              Constraint: LDQ .ge. N.
C
C  8: ERR    - REAL.                                              Output
C
C              On exit: a measure of the accuracy  of the eigenvalues of
C              T compared to the requested eigenvalues,  ERR returns the
C              value err defined in Section 7.
C
C  9: WORK(*) - REAL array.                             Workspace/Output
C
C              Note:  the dimension of the array  WORK  must be at least
C              N*(N + 2).
C
C              On exit: the first element of WORK returns the number of
C              iterations.
C
C 10: IFAIL  - INTEGER.                                     Input/Output
C
C              On entry:  IFAIL must be set to 0, -1 or 1. For users not
C              familiar  with this  parameter (described in Chapter P01)
C              the recommended value is 0.
C
C              On exit:  IFAIL = 0  unless  the routine detects an error
C              (see Section 6).
C
C     5. Error Indicators and Warnings
C     ================================
C
C     Errors detected by the routine:
C
C     If on entry  IFAIL = 0 or -1,  explanatory error messages  are
C     output on the current error message unit (as defined by X04AAF).
C
C     IFAIL = -1
C
C     One or more of the following conditions holds:
C
C        N   .lt. 1
C        LDQ .lt. N
C
C     IFAIL .eq. 1
C
C     The iterative process  has failed to converge in  50*N iterations.
C     The  current approximation  to  the  parameters  of  the  Toeplitz
C     matrix are returned in T, and ERR returns the value err defined in
C     Section 7.
C
C     6. Accuracy
C     ===========
C
C  The  (computed)  eigenvalues   of  the   Toeplitz  matrix  T  satisfy
C
C     err .le. N*eps,
C
C     where err is defined as
C
C     err = max( abs( s(i) - d(i) ) )/( 1 + max( abs( d(i) ) ) ),
C
C  s(i) is the i(th) eigenvalue of  T and  eps is the machine precision.
C
C     7. Further Comments
C     ===================
C
C     None.
C
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='Y90RQF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ERR, T0
      INTEGER           IFAIL, JOB, LDQ, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(N), Q(LDQ,N), T(*), WORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  DMAX, TOL
      INTEGER           IERR, ITER, ITMAX, J
C     .. Local Arrays ..
      CHARACTER*48      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06FBF, M01CAF, P01ABY, Y90RQX, Y90RQY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Check the input parameters.
C
C-----------------------------------------------------------------------
      IERR = 0
      IF (N.LT.1) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (LDQ.LT.N) CALL P01ABY(LDQ,'LDQ',IFAIL,IERR,SRNAME)
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IFAIL = P01ABF(IFAIL,-1,SRNAME,1,REC)
         RETURN
      END IF
C-----------------------------------------------------------------------
C
C     Set the  convergence tolerance  and  maximum number of iterations.
C
C-----------------------------------------------------------------------
      IF (JOB.GE.3) THEN
         TOL = N*X02AJF()
         ITMAX = 50*N
      ELSE
         TOL = N*SQRT(X02AJF())
         ITMAX = 5*N
      END IF
C-----------------------------------------------------------------------
C
C     Sort the eigenvalues into ascending order and find the eigenvalue
C     of largest absolute value.
C
C-----------------------------------------------------------------------
      IERR = 1
      CALL M01CAF(D,1,N,'Ascending',IERR)
      DMAX = MAX(ABS(D(1)),ABS(D(N)))
C-----------------------------------------------------------------------
C
C     Initialize T and Q
C
C-----------------------------------------------------------------------
      IF ((JOB.LE.1) .OR. (JOB.EQ.4)) THEN
C
C        Initialize T and S (first n elements of WORK),  in case of
C        failure in  either  F04ARF or F02ABF  (which are called from
C        Y90RQX),  in which case we simply run out of iterations.
C
         CALL F06FBF(N,ZERO,T,1)
         CALL F06FBF(N,ZERO,WORK,1)
C
C        Form the initial  Q  matrix.
C
         CALL Y90RQY(N,Q,LDQ)
C
      END IF
C-----------------------------------------------------------------------
C
C     Iterate ...
C
C-----------------------------------------------------------------------
      ITER = 0
   20 CONTINUE
C-----------------------------------------------------------------------
C
C     repeat
C
C-----------------------------------------------------------------------
      ITER = ITER + 1
      CALL Y90RQX(N,Q,LDQ,D,T0,T,WORK,WORK(N+1),N,WORK(N**2+N+1),IERR)
C-----------------------------------------------------------------------
C
C        Test for convergence.
C
C-----------------------------------------------------------------------
      ERR = ZERO
      DO 40 J = 1, N
         ERR = MAX(ERR,ABS(D(J)-WORK(J)))
   40 CONTINUE
      ERR = ERR/(1+DMAX)
      IF ((ERR.GT.TOL) .AND. (ITER.LT.ITMAX)) GO TO 20
C-----------------------------------------------------------------------
C
C     until( ( ERR.LE.TOL ).OR. ( ITER.EQ.50*N ) )
C
C-----------------------------------------------------------------------
      WORK(1) = ITER
      IF ((ITER.EQ.ITMAX) .AND. (JOB.GE.3)) THEN
         WRITE (REC,FMT=99998)
         IFAIL = P01ABF(IFAIL,1,SRNAME,1,REC)
      ELSE
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90RQF.
C
C-----------------------------------------------------------------------
      RETURN
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
99998 FORMAT ('    The iterative process has failed to converge')
      END
