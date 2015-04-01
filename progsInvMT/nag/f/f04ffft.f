      SUBROUTINE F04FFF(N,T,B,X,WANTP,P,WORK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     -- Written on 1-December-1990.
C     This version dated 29-December-1990.
C     Sven Hammarling, Nag Ltd.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04FFF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, N
      LOGICAL           WANTP
C     .. Array Arguments ..
      DOUBLE PRECISION  B(*), P(*), T(0:*), WORK(*), X(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, BETA
      INTEGER           I, IERR
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
      IF (N.EQ.0) THEN
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
         RETURN
      END IF
C
C     Start the recursion.
C
      BETA = T(0)
      X(1) = B(1)/BETA
      DO 20 I = 2, N
C
C        work(1), ..., work(i-2)   and  work(n), ..., work(n+i-3)   both
C        contain  the  solution  to the  Yule-Walker equations  of order
C        (i-2).  (The second copy is in reverse order.)  We  first solve
C        the Yule-Walker equations of order (i-1).
C
         ALPHA = -(T(I-1)+DDOT(I-2,T(1),1,WORK(N),1))/BETA
         IF (WANTP) P(I-1) = ALPHA
         IF (ABS(ALPHA).GE.ONE) THEN
            WRITE (REC,FMT=99998) I, ALPHA
            IFAIL = P01ABF(IFAIL,I,SRNAME,2,REC)
            RETURN
         END IF
         CALL DAXPY(I-2,ALPHA,WORK(N),1,WORK,1)
         WORK(I-1) = ALPHA
         BETA = (ONE-ALPHA)*(ONE+ALPHA)*BETA
         CALL DCOPY(I-1,WORK,-1,WORK(N),1)
C
C        Now form  the  solution for the equations of  order i  with the
C        general right hand side.
C
         X(I) = (B(I)-DDOT(I-1,T(1),1,X,-1))/BETA
         CALL DAXPY(I-1,X(I),WORK(N),1,X,1)
   20 CONTINUE
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C
C     End of F04FFF.
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
99998 FORMAT ('    Principal minor ',I8,' is not positive-definite',
     *       /'    Value of the reflection coefficient is ',1P,D10.2)
      END
