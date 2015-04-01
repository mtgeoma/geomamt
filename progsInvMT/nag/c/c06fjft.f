      SUBROUTINE C06FJF(NDIM,ND,N,X,Y,WORK,LWORK,IFAIL)
C     MARK 11 RELEASE. NAG COPYRIGHT 1984.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     MULTI-DIMENSIONAL DISCRETE FOURIER TRANSFORM OF A
C     MULTI-DIMENSIONAL SEQUENCE OF COMPLEX DATA VALUES
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06FJF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LWORK, N, NDIM
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(LWORK), X(N), Y(N)
      INTEGER           ND(NDIM)
C     .. Local Scalars ..
      INTEGER           IFAIL1, L
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C06FFF
C     .. Executable Statements ..
      IF (NDIM.LT.1) GO TO 40
      DO 20 L = 1, NDIM
         IFAIL1 = 1
         CALL C06FFF(NDIM,L,ND,N,X,Y,WORK,LWORK,IFAIL1)
         IF (IFAIL1.NE.0) GO TO 60
   20 CONTINUE
      IFAIL = 0
      RETURN
   40 IFAIL1 = 1
   60 IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,0,P01REC)
      RETURN
      END
