      SUBROUTINE F07AUZ(N,SA,CX,INCX)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             ZDRSCL(N,SA,CX,INCX)
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
C     Courant Institute, NAG Ltd., and Univ. of California Berkeley
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SA
      INTEGER           INCX, N
C     .. Array Arguments ..
      COMPLEX*16        CX(*)
C     .. Local Scalars ..
      INTEGER           I, IX, M, MP1
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
C
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 40
C
C     Code for increment not equal to 1
C
      IX = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      DO 20 I = 1, N
         CX(IX) = CX(IX)/SA
         IX = IX + INCX
   20 CONTINUE
      RETURN
C
C     Code for increment equal to 1
C
C
C     Clean-up loop
C
   40 CONTINUE
      M = MOD(N,5)
      IF (M.EQ.0) GO TO 80
      DO 60 I = 1, M
         CX(I) = CX(I)/SA
   60 CONTINUE
      IF (N.LT.5) RETURN
   80 CONTINUE
      MP1 = M + 1
      DO 100 I = MP1, N, 5
         CX(I) = CX(I)/SA
         CX(I+1) = CX(I+1)/SA
         CX(I+2) = CX(I+2)/SA
         CX(I+3) = CX(I+3)/SA
         CX(I+4) = CX(I+4)/SA
  100 CONTINUE
      RETURN
      END
