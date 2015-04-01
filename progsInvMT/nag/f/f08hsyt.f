      SUBROUTINE F08HSY(N,X,INCX,Y,INCY,C,S,INCC)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZLARTV(N,X,INCX,Y,INCY,C,S,INCC)
C
C  Purpose
C  =======
C
C  ZLARTV applies a vector of complex plane rotations with real cosines
C  to elements of the complex vectors x and y. For i = 1,2,...,n
C
C     ( x(i) ) := (        c(i)   s(i) ) ( x(i) )
C     ( y(i) )    ( -conjg(s(i))  c(i) ) ( y(i) )
C
C  Arguments
C  =========
C
C  N       (input) INTEGER
C          The number of plane rotations to be applied.
C
C  X       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCX)
C          The vector x.
C
C  INCX    (input) INTEGER
C          The increment between elements of X. INCX > 0.
C
C  Y       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCY)
C          The vector y.
C
C  INCY    (input) INTEGER
C          The increment between elements of Y. INCY > 0.
C
C  C       (input) DOUBLE PRECISION array, dimension (1+(N-1)*INCC)
C          The cosines of the plane rotations.
C
C  S       (input) COMPLEX*16 array, dimension (1+(N-1)*INCC)
C          The sines of the plane rotations.
C
C  INCC    (input) INTEGER
C          The increment between elements of C and S. INCC > 0.
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Scalar Arguments ..
      INTEGER           INCC, INCX, INCY, N
C     .. Array Arguments ..
      COMPLEX*16        S(*), X(*), Y(*)
      DOUBLE PRECISION  C(*)
C     .. Local Scalars ..
      COMPLEX*16        XI, YI
      INTEGER           I, IC, IX, IY
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG
C     .. Executable Statements ..
C
      IX = 1
      IY = 1
      IC = 1
      DO 20 I = 1, N
         XI = X(IX)
         YI = Y(IY)
         X(IX) = C(IC)*XI + S(IC)*YI
         Y(IY) = C(IC)*YI - DCONJG(S(IC))*XI
         IX = IX + INCX
         IY = IY + INCY
         IC = IC + INCC
   20 CONTINUE
      RETURN
C
C     End of F08HSY (ZLARTV)
C
      END
