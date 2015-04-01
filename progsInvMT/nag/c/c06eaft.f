      SUBROUTINE C06EAF(X,PTS,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     REAL FOURIER TRANSFORM
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06EAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X(PTS)
C     .. Local Scalars ..
      DOUBLE PRECISION  SQPTS
      INTEGER           IERROR, IPTS, PMAX, PSYM, TWOGRP
C     .. Local Arrays ..
      INTEGER           FACTOR(21), SYM(21), UNSYM(21)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C06EAW, C06EAY, C06EAZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Data statements ..
      DATA              PMAX/19/
      DATA              TWOGRP/8/
C     .. Executable Statements ..
      IF (PTS.LE.1) GO TO 40
      IERROR = 0
      CALL C06EAZ(PTS,PMAX,TWOGRP,FACTOR,SYM,PSYM,UNSYM,IERROR)
      IF (IERROR.NE.0) GO TO 60
      CALL C06EAY(X,PTS,SYM,PSYM,UNSYM)
      CALL C06EAW(X,PTS,FACTOR)
      SQPTS = SQRT(DBLE(PTS))
      DO 20 IPTS = 1, PTS
         X(IPTS) = X(IPTS)/SQPTS
   20 CONTINUE
      IFAIL = 0
      GO TO 80
C
   40 IERROR = 3
   60 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
   80 RETURN
      END
