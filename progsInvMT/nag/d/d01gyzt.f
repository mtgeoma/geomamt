      INTEGER FUNCTION D01GYZ(N)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RETURNS 0 IF N IS PRIME
C     RETURNS 1 IF N IS NOT PRIME
C     .. Scalar Arguments ..
      INTEGER                 N
C     .. Local Scalars ..
      INTEGER                 ND, NINC
C     .. Intrinsic Functions ..
      INTRINSIC               MOD
C     .. Executable Statements ..
      D01GYZ = 1
      IF (N.EQ.2 .OR. N.EQ.3 .OR. N.EQ.5) GO TO 40
      IF (MOD(N,2).EQ.0) RETURN
      IF (MOD(N,3).EQ.0) RETURN
      ND = 5
      NINC = 2
   20 IF (MOD(N,ND).EQ.0) RETURN
      ND = ND + NINC
      NINC = 6 - NINC
      IF (ND.LE.N/ND) GO TO 20
   40 D01GYZ = 0
      RETURN
      END
