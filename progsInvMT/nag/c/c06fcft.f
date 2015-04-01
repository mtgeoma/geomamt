      SUBROUTINE C06FCF(X,Y,PTS,WORK,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     COMPLEX FOURIER TRANSFORM
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06FCF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(PTS), X(PTS), Y(PTS)
C     .. Local Scalars ..
      DOUBLE PRECISION  SQPTS
      INTEGER           IERROR, IPTS, PMAX, TWOGRP
C     .. Local Arrays ..
      INTEGER           RFACT(21), TFACT(21)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C06ECW, C06FAY, C06FAZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Data statements ..
      DATA              PMAX/19/
      DATA              TWOGRP/8/
C     .. Executable Statements ..
      IF (PTS.LE.1) GO TO 40
      CALL C06FAZ(PTS,PMAX,TWOGRP,TFACT,RFACT,IERROR)
      IF (IERROR.NE.0) GO TO 60
      CALL C06ECW(X,Y,PTS,TFACT)
      CALL C06FAY(X,PTS,RFACT,WORK)
      CALL C06FAY(Y,PTS,RFACT,WORK)
      SQPTS = SQRT(DBLE(PTS))
      DO 20 IPTS = 1, PTS
         X(IPTS) = X(IPTS)/SQPTS
         Y(IPTS) = Y(IPTS)/SQPTS
   20 CONTINUE
      IFAIL = 0
      GO TO 80
C
   40 IERROR = 3
   60 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
   80 RETURN
      END
