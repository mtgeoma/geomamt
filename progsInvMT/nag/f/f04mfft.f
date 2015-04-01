      SUBROUTINE F04MFF(N,T,B,X,P,WORK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     -- Written on 9-February-1990.
C     This version dated 8-January-1991.
C     Sven Hammarling, Nag Ltd.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04MFF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  P
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  B(*), T(0:*), WORK(*), X(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  BETA
      INTEGER           IERR
C     .. Local Arrays ..
      CHARACTER*53      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      INTEGER           P01ABF
      EXTERNAL          DDOT, P01ABF
C     .. External Subroutines ..
      EXTERNAL          P01ABX, P01ABY, DAXPY, DCOPY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
C
C
C     Check the input parameters.
C
      IERR = 0
      IF (N.LT.0) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (T(0).LE.ZERO) CALL P01ABX(T(0),'T(0)',IFAIL,IERR,SRNAME)
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IFAIL = P01ABF(IFAIL,-1,SRNAME,1,REC)
         RETURN
      END IF
C
C     Solve the equations.  On exit the element work(2*n-1) stores beta,
C     work(1), ..., work(n-1)   and   work(n), ..., work(2*n-2)  contain
C     the solution  to the  Yule-Walker equations  of  order (n-1). (The
C     second copy is in reverse order.)
C
      IF (N.EQ.1) THEN
         X(1) = B(1)/T(0)
         WORK(1) = T(0)
      ELSE IF (N.GT.1) THEN
C
C        First solve the Yule-Walker equations.
C
         BETA = WORK(2*N-3)
         P = -(T(N-1)+DDOT(N-2,T(1),1,WORK(N-1),1))/BETA
         IF (ABS(P).GE.ONE) THEN
            WRITE (REC,FMT=99998) P
            IFAIL = P01ABF(IFAIL,1,SRNAME,2,REC)
            RETURN
         END IF
         CALL DAXPY(N-2,P,WORK(N-1),1,WORK,1)
         WORK(N-1) = P
         BETA = (ONE-P)*(ONE+P)*BETA
         CALL DCOPY(N-1,WORK,-1,WORK(N),1)
C
C        Now solve for the general right-hand side.
C
         X(N) = (B(N)-DDOT(N-1,T(1),1,X,-1))/BETA
         CALL DAXPY(N-1,X(N),WORK(N),1,X,1)
         WORK(2*N-1) = BETA
      END IF
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C
C     End of F04MFF.
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
99998 FORMAT ('    The Toeplitz Matrix is not positive-definite',/'   ',
     *       ' Value of the reflection coefficient is ',1P,D10.2)
      END
