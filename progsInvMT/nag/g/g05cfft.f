      SUBROUTINE G05CFF(IA,NI,XA,NX,IFAIL)
C     MARK 14 RE-ISSUE. NAG COPYRIGHT 1989.
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     THIS SAVES THE GENERATOR STATE IN IA AND XA.
C
C     This revised version of G05CFF has been introduced for
C     compatibility with the new routine G05FDF, introduced at Mark 14.
C     G05CFZ now saves three values in XA(2), XA(3) and XA(4) for use
C     by G05FDF, G05DDF and G05DGF respectively, and XA(1) is used for
C     check-summing.
C
C     Jeremy Du Croz, NAG Ltd, June 1989.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G05CFF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, NI, NX
C     .. Array Arguments ..
      DOUBLE PRECISION  XA(NX)
      INTEGER           IA(NI)
C     .. Local Scalars ..
      DOUBLE PRECISION  FACT, X, ZERO
      INTEGER           I, IERR, J
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G05CFZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Data statements ..
      DATA              FACT/4096.0D0/, ZERO/0.0D0/
C     .. Executable Statements ..
      IERR = 1
      IF (NI.LT.9) GO TO 60
      IERR = 2
      IF (NX.LT.4) GO TO 60
      CALL G05CFZ(IA(2),XA(2))
      J = 0
      DO 20 I = 2, 9
         J = J + IA(I)
   20 CONTINUE
      IA(1) = J
      X = DBLE(J)/FACT
      DO 40 I = 2, 4
         X = X + XA(I)
   40 CONTINUE
      XA(1) = X
      IFAIL = 0
      RETURN
   60 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
      END
