      SUBROUTINE G05EJF(IA,N,IZ,M,IFAIL)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G05EJF TAKES A RANDOM SAMPLE OF SIZE M FROM IA (OF
C     SIZE N) AND PUTS IT INTO IZ.
C
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G05EJF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M, N
C     .. Array Arguments ..
      INTEGER           IA(N), IZ(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  ONE, P, Q
      INTEGER           I, IERR, J
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G05CAF
      INTEGER           P01ABF
      EXTERNAL          G05CAF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Data statements ..
      DATA              ONE/1.0D0/
C     .. Executable Statements ..
      IERR = 1
      IF (N.LT.1) GO TO 60
      IERR = 2
      IF ((M.LT.1) .OR. (M.GT.N)) GO TO 60
      P = DBLE(M)
      Q = DBLE(N)
      J = 1
      DO 40 I = 1, N
         IF (Q*G05CAF(0.0D0).GT.P) GO TO 20
         IZ(J) = IA(I)
         J = J + 1
         P = P - ONE
   20    Q = Q - ONE
   40 CONTINUE
      IFAIL = 0
      RETURN
   60 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
      END
