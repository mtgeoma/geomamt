      SUBROUTINE F04FEF(N,T,X,WANTP,P,WANTV,V,VLAST,WORK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     -- Written on 9-February-1990.
C     This version dated 29-December-1990.
C     Sven Hammarling, Nag Ltd.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04FEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  VLAST
      INTEGER           IFAIL, N
      LOGICAL           WANTP, WANTV
C     .. Array Arguments ..
      DOUBLE PRECISION  P(*), T(0:N), V(*), WORK(*), X(*)
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
C     Start the recursion.
C
      BETA = T(0)
      DO 20 I = 1, N
C
C        Make  a  copy  of the  previous solution (in reverse order)  to
C        simplify the code and to allow the call to DAXPY.
C
         CALL DCOPY(I-1,X,-1,WORK,1)
         ALPHA = -(T(I)+DDOT(I-1,T(1),1,WORK,1))/BETA
         CALL DAXPY(I-1,ALPHA,WORK,1,X,1)
         X(I) = ALPHA
         BETA = (ONE-ALPHA)*(ONE+ALPHA)*BETA
         IF (WANTP) P(I) = ALPHA
         IF (WANTV) V(I) = BETA/T(0)
         IF (ABS(ALPHA).GE.ONE) THEN
            VLAST = BETA/T(0)
            WRITE (REC,FMT=99998) I + 1, ALPHA
            IFAIL = P01ABF(IFAIL,I,SRNAME,2,REC)
            RETURN
         END IF
   20 CONTINUE
      VLAST = BETA/T(0)
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C
C     End of F04FEF.
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
99998 FORMAT ('    Principal minor ',I8,' is not positive-definite',
     *       /'    Value of the reflection coefficient is ',1P,D10.2)
      END
