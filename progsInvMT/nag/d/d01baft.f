      DOUBLE PRECISION FUNCTION D01BAF(WTFUN,A,B,NPTS,FUN,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 8A REVISED. IER-250 (JULY 1980).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     EVALUATES INTEGRAL OF FUNCTION FUN BY NPTS GAUSS FORMULA OF
C     TYPE WTFUN
C     IFAIL = 1 - THE NPTS RULE IS NOT AMONG THOSE STORED
C     ( ANSWER EVALUATED FOR LARGEST VALID NPTS LESS THAN REQUESTED
C     VALUE)
C     IFAIL = 2 - VALUES OF A OR B INVALID
C     ( ANSWER RETURNED AS ZERO)
C
C     WTFUN
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='D01BAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, B
      INTEGER                          IFAIL, NPTS
C     .. Function Arguments ..
      DOUBLE PRECISION                 FUN
      EXTERNAL                         FUN
C     .. Subroutine Arguments ..
      EXTERNAL                         WTFUN
C     .. Local Scalars ..
      DOUBLE PRECISION                 AA, SUM, WW
      INTEGER                          I, IERR, NPTSA
C     .. Local Arrays ..
      DOUBLE PRECISION                 ABSCIS(64), WEIGHT(64)
      CHARACTER*1                      P01REC(1)
C     .. External Functions ..
      INTEGER                          P01ABF
      EXTERNAL                         P01ABF
C     .. External Subroutines ..
      EXTERNAL                         D01BBF
C     .. Intrinsic Functions ..
      INTRINSIC                        MIN
C     .. Executable Statements ..
      SUM = 0.0D0
      NPTSA = MIN(NPTS,64)
      IERR = 1
      CALL D01BBF(WTFUN,A,B,1,NPTSA,WEIGHT,ABSCIS,IERR)
      IF (IERR.GT.1) GO TO 40
      IF (NPTS.NE.NPTSA) IERR = 1
      DO 20 I = 1, NPTSA
         WW = WEIGHT(I)
         IF (WW.EQ.0.0D0) GO TO 20
         AA = ABSCIS(I)
         SUM = SUM + WW*FUN(AA)
   20 CONTINUE
      IF (IERR.EQ.0) GO TO 60
   40 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      GO TO 80
   60 IFAIL = 0
   80 D01BAF = SUM
      RETURN
      END
