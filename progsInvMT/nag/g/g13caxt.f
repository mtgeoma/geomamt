      SUBROUTINE G13CAX(XG,NX,M,PI)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G13CAX  APPLIES THE SPLIT COSINE BELL TAPER
C
C     XG    - DATA ARRAY
C     NX    - NUMBER OF DATA POINTS
C     M     - NUMBER OF POINTS TAPERED AT EACH END
C     PI    - CONSTANT PI
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  PI
      INTEGER           M, NX
C     .. Array Arguments ..
      DOUBLE PRECISION  XG(NX)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, S1
      INTEGER           I, NI
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE
C     .. Executable Statements ..
      IF (M.LE.0) GO TO 40
      A = DBLE(M)
      NI = NX + 1
      DO 20 I = 1, M
         S1 = COS(PI*(DBLE(I)-0.5D0)/A)
         S1 = (1.0D0-S1)/2.0D0
         NI = NI - 1
         XG(I) = XG(I)*S1
         XG(NI) = XG(NI)*S1
   20 CONTINUE
   40 CONTINUE
      RETURN
      END
