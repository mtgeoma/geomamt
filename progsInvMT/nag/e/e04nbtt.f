      SUBROUTINE E04NBT(MODE,NROWT,N,T,Y)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-587 (MAR 1988).
C
C     ******************************************************************
C     E04NBT  solves equations involving a reverse-triangular matrix  T
C     and a right-hand-side vector  y,  returning the solution in  y.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original Fortran 77 version written February-1985.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           MODE, N, NROWT
C     .. Array Arguments ..
      DOUBLE PRECISION  T(NROWT,*), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  YJ
      INTEGER           J, JJ, L, N1
C     .. External Subroutines ..
      EXTERNAL          DAXPY
C     .. Executable Statements ..
C
      N1 = N + 1
      IF (MODE.EQ.1) THEN
C
C        Mode = 1  ---  Solve  T * y(new) = y(old).
C
         DO 20 J = 1, N
            JJ = N1 - J
            YJ = Y(J)/T(J,JJ)
            Y(J) = YJ
            L = JJ - 1
            IF (L.GT.0 .AND. YJ.NE.ZERO) CALL DAXPY(L,(-YJ),T(J+1,JJ),1,
     *          Y(J+1),1)
   20    CONTINUE
      ELSE
C
C        Mode = 2  ---  Solve  T' y(new) = y(old).
C
         DO 40 J = 1, N
            JJ = N1 - J
            YJ = Y(J)/T(JJ,J)
            Y(J) = YJ
            L = JJ - 1
            IF (L.GT.0 .AND. YJ.NE.ZERO) CALL DAXPY(L,(-YJ),T(JJ,J+1),
     *          NROWT,Y(J+1),1)
   40    CONTINUE
      END IF
C
C     Reverse the solution vector.
C
      IF (N.GT.1) THEN
         L = N/2
         DO 60 J = 1, L
            JJ = N1 - J
            YJ = Y(J)
            Y(J) = Y(JJ)
            Y(JJ) = YJ
   60    CONTINUE
      END IF
C
      RETURN
C
C     End of  E04NBT. (CMTSOL)
C
      END
