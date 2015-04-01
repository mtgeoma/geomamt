      SUBROUTINE G04CAW(N,NLEVEL,B,T,SEM,ESS,R,IFAC)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Computes treatment sums of squares, treatment means and s.e.s
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ESS, SEM
      INTEGER           N, NLEVEL
C     .. Array Arguments ..
      DOUBLE PRECISION  B(NLEVEL), R(N), T(NLEVEL)
      INTEGER           IFAC(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  SCALE, TEMP
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      SCALE = N/NLEVEL
      SCALE = 1.0D0/SCALE
      SEM = SQRT(2.0D0*SCALE)
      ESS = 0.0D0
      DO 20 I = 1, NLEVEL
         TEMP = B(I)
         ESS = ESS + TEMP*TEMP
         B(I) = TEMP*SCALE
         T(I) = T(I)*SCALE
   20 CONTINUE
      ESS = ESS*SCALE
      DO 40 I = 1, N
         R(I) = R(I) - B(IFAC(I))
   40 CONTINUE
      RETURN
      END
