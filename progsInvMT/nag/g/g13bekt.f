      SUBROUTINE G13BEK(A,IA,B,IB,APB,IAPB,NA,NB,NAPB)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12 REVISED. IER-522 (AUG 1986).
C     MARK 14B REVISED. IER-842 (MAR 1990).
C
C     ROUTINE G13BEK TAKES NA VALUES OF ARRAY A AND NB VALUES
C     OF ARRAY B AND COMBINES THEM IN THAT ORDER IN A NEW ARRAY APB
C
C     .. Scalar Arguments ..
      INTEGER           IA, IAPB, IB, NA, NAPB, NB
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA), APB(IAPB), B(IB)
C     .. Local Scalars ..
      INTEGER           I, J, K, NPB
C     .. Executable Statements ..
      NAPB = NA + NB
      IF (NB.LE.0) GO TO 60
      DO 20 I = 1, NB
         APB(I) = B(I)
   20 CONTINUE
      IF (NA.LE.0) GO TO 100
      NPB = NB + 1
      DO 40 I = 1, NB
         J = NPB - I
         K = NA + J
         APB(K) = APB(J)
   40 CONTINUE
   60 IF (NA.LE.0) GO TO 100
      DO 80 I = 1, NA
         APB(I) = A(I)
   80 CONTINUE
  100 RETURN
      END
