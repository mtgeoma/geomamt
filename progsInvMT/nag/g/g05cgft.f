      SUBROUTINE G05CGF(IA,NI,XA,NX,IFAIL)
C     MARK 14 RE-ISSUE. NAG COPYRIGHT 1989.
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     THIS RESTORES THE GENERATOR STATE FROM IA AND XA.
C
C     This revised version of G05CGF has been introduced for
C     compatibility with the new routine G05FDF, introduced at Mark 14.
C     G05CGZ now restores three values from XA(2), XA(3) and XA(4) for
C     use by G05FDF, G05DDF and G05DGF respectively, and XA(1) is used
C     for check-summing.
C
C     Jeremy Du Croz, NAG Ltd, June 1989.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G05CGF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, NI, NX
C     .. Array Arguments ..
      DOUBLE PRECISION  XA(NX)
      INTEGER           IA(NI)
C     .. Local Scalars ..
      DOUBLE PRECISION  FACT, X, ZERO
      INTEGER           I, IERR, J
C     .. Local Arrays ..
      CHARACTER         P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G05CGZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Data statements ..
      DATA              FACT/4096.0D0/, ZERO/0.0D0/
C     .. Executable Statements ..
      IERR = 1
      IF ((NI.LT.9) .OR. (NX.LT.4)) GO TO 60
      J = 0
      DO 20 I = 2, 9
         J = J + IA(I)
   20 CONTINUE
      IERR = 2
      IF (J.NE.IA(1)) GO TO 60
      X = DBLE(J)/FACT
      DO 40 I = 2, 4
         X = X + XA(I)
   40 CONTINUE
      IF (ABS(XA(1)-X).GE.4.0D0*X02AJF()*ABS(XA(1))) GO TO 60
      CALL G05CGZ(IA(2),XA(2),I)
      IERR = 3
      IF (I.NE.0) GO TO 60
      IFAIL = 0
      RETURN
   60 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
      END
