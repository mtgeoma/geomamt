      DOUBLE PRECISION FUNCTION D01FDW(F,Y,NDIM,REGION)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     TRANSFORMATION ROUTINE
C     REGION
C     .. Scalar Arguments ..
      INTEGER                          NDIM
C     .. Array Arguments ..
      DOUBLE PRECISION                 Y(NDIM)
C     .. Function Arguments ..
      DOUBLE PRECISION                 F
      EXTERNAL                         F
C     .. Subroutine Arguments ..
      EXTERNAL                         REGION
C     .. Local Scalars ..
      DOUBLE PRECISION                 C, D, FACT
      INTEGER                          J
C     .. Local Arrays ..
      DOUBLE PRECISION                 YP(30)
C     .. Executable Statements ..
      FACT = 1.0D0
      DO 20 J = 1, NDIM
         CALL REGION(NDIM,YP,J,C,D)
         FACT = FACT*(D-C)*0.5D0
         YP(J) = (D+C+(D-C)*Y(J))*0.5D0
   20 CONTINUE
      D01FDW = F(NDIM,YP)*FACT
      RETURN
      END
