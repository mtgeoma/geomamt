      SUBROUTINE F04MEF(N,T,X,V,WORK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     -- Written on 9-February-1990.
C     This version dated 8-January-1991.
C     Sven Hammarling, Nag Ltd.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04MEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  V
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  T(0:N), WORK(*), X(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA
      INTEGER           IERR
C     .. Local Arrays ..
      CHARACTER*58      REC(2)
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
      IF ((N.GT.1) .AND. (ABS(X(N-1)).GE.ONE)) CALL P01ABX(X(N-1),
     *    'X(N-1)',IFAIL,IERR,SRNAME)
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IFAIL = P01ABF(IFAIL,-1,SRNAME,1,REC)
         RETURN
      END IF
C
      IF (N.GT.0) THEN
         IF (N.EQ.1) V = ONE
C
C        We  make  a  copy  (in reverse order)  of the  solution  of the
C        equations of  order (n-1).  This simplifies the code and allows
C        us to call DAXPY.
C
         CALL DCOPY(N-1,X,-1,WORK,1)
         ALPHA = -(T(N)+DDOT(N-1,T(1),1,WORK,1))/(V*T(0))
         CALL DAXPY(N-1,ALPHA,WORK,1,X,1)
         X(N) = ALPHA
         V = (ONE-ALPHA)*(ONE+ALPHA)*V
         IF (ABS(ALPHA).GE.ONE) THEN
            WRITE (REC,FMT=99998) N + 1, ALPHA
            IFAIL = P01ABF(IFAIL,1,SRNAME,2,REC)
            RETURN
         END IF
      END IF
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C
C     End of F04MEF.
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
99998 FORMAT ('    Matrix of order',I8,
     *       ' would not be positive-definite',/'    Value of the refl',
     *       'ection coefficient is ',1P,D10.2)
      END
