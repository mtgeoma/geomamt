      SUBROUTINE F08HEZ(N,X,INCX,Y,INCY,C,INCC)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DLARGV(N,X,INCX,Y,INCY,C,INCC)
C
C  Purpose
C  =======
C
C  DLARGV generates a vector of real plane rotations, determined by
C  elements of the real vectors x and y. For i = 1,2,...,n
C
C     (  c(i)  s(i) ) ( x(i) ) = ( a(i) )
C     ( -s(i)  c(i) ) ( y(i) ) = (   0  )
C
C  Arguments
C  =========
C
C  N       (input) INTEGER
C          The number of plane rotations to be generated.
C
C  X       (input/output) DOUBLE PRECISION array,
C                         dimension (1+(N-1)*INCX)
C          On entry, the vector x.
C          On exit, x(i) is overwritten by a(i), for i = 1,...,n.
C
C  INCX    (input) INTEGER
C          The increment between elements of X. INCX > 0.
C
C  Y       (input/output) DOUBLE PRECISION array,
C                         dimension (1+(N-1)*INCY)
C          On entry, the vector y.
C          On exit, the sines of the plane rotations.
C
C  INCY    (input) INTEGER
C          The increment between elements of Y. INCY > 0.
C
C  C       (output) DOUBLE PRECISION array, dimension (1+(N-1)*INCC)
C          The cosines of the plane rotations.
C
C  INCC    (input) INTEGER
C          The increment between elements of C. INCC > 0.
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INCC, INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION  C(*), X(*), Y(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  TT, W, XI, YI
      INTEGER           I, IC, IX, IY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Executable Statements ..
C
      IX = 1
      IY = 1
      IC = 1
      DO 20 I = 1, N
         XI = X(IX)
         YI = Y(IY)
         IF (XI.EQ.ZERO) THEN
            C(IC) = ZERO
            Y(IY) = ONE
            X(IX) = YI
         ELSE
            W = MAX(ABS(XI),ABS(YI))
            XI = XI/W
            YI = YI/W
            TT = SQRT(XI*XI+YI*YI)
            C(IC) = XI/TT
            Y(IY) = YI/TT
            X(IX) = W*TT
         END IF
         IX = IX + INCX
         IY = IY + INCY
         IC = IC + INCC
   20 CONTINUE
      RETURN
C
C     End of F08HEZ (DLARGV)
C
      END
