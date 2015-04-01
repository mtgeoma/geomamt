      SUBROUTINE C06FFF(NDIM,L,ND,N,X,Y,WORK,LWORK,IFAIL)
C     MARK 11 RELEASE. NAG COPYRIGHT 1984.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     DISCRETE FOURIER TRANSFORM OF ONE VARIABLE IN A
C     MULTI-VARIABLE SEQUENCE OF COMPLEX DATA VALUES
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06FFF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, L, LWORK, N, NDIM
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(LWORK), X(N), Y(N)
      INTEGER           ND(NDIM)
C     .. Local Scalars ..
      INTEGER           I, IFAIL1, NI, NK, NL
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C06FFZ
C     .. Executable Statements ..
      IF (NDIM.LT.1) GO TO 60
      IF (L.LT.1 .OR. L.GT.NDIM) GO TO 100
      NK = 1
      NI = 1
      DO 20 I = 1, NDIM
         IF (ND(I).LT.1) GO TO 120
         IF (I.GT.L) NK = NK*ND(I)
         IF (I.LT.L) NI = NI*ND(I)
   20 CONTINUE
      NL = ND(L)
      IF (NI*NL*NK.NE.N) GO TO 80
      IF (NL.EQ.1) GO TO 40
      IF (LWORK.LT.3*NL) GO TO 160
      CALL C06FFZ(X,Y,NI,NL,NK,WORK,WORK(NL+1),WORK(2*NL+1),IFAIL1)
      IF (IFAIL1.NE.0) GO TO 140
   40 IFAIL = 0
      RETURN
   60 IFAIL1 = 1
      GO TO 180
   80 IFAIL1 = 2
      GO TO 180
  100 IFAIL1 = 3
      GO TO 180
  120 IFAIL1 = 3 + 10*I
      GO TO 180
  140 IFAIL1 = IFAIL1 + 10*L
      GO TO 180
  160 IFAIL1 = 4 + 10*L
  180 IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,0,P01REC)
      RETURN
      END
