      SUBROUTINE G07EBX(N,X,M,Y,D,U,TOL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     This routine computes the Mann-Whitney U test statistic
C     for the two samples X(i) + D and Y(j).
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  D, TOL, U
      INTEGER           M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N), Y(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP, TU, XI, XTNI, XTPI
      INTEGER           I, IU, J, JLE, K
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
      JLE = 0
      IU = 0
      TU = 0.0D0
      U = 0.0D0
C
C     Main loop
C
      DO 40 I = 1, N
         XI = X(I) + D
         XTNI = XI - TOL
         XTPI = XI + TOL
         J = JLE
         K = 0
   20    CONTINUE
         TEMP = Y(J+1)
         IF (XTNI.GE.TEMP) THEN
            JLE = JLE + 1
            J = J + 1
            IF (JLE.GE.M) THEN
C
C              X(I) is greater than (or equal to) all Y(J). Therefore
C              X(I+1),...,X(N) are all greater than all Y(J).
C
               IU = IU + (N-I+1)*M
               GO TO 60
            END IF
            GO TO 20
         ELSE IF (XTPI.GT.TEMP) THEN
            K = K + 1
            J = J + 1
            IF (J.LT.M) GO TO 20
         END IF
         IU = IU + JLE
         TU = TU + DBLE(K)/2.0D0
   40 CONTINUE
C
   60 U = DBLE(IU) + TU
      RETURN
      END
