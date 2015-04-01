      SUBROUTINE F08HEX(N,X,Y,Z,INCX,C,S,INCC)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DLAR2V(N,X,Y,Z,INCX,C,S,INCC)
C
C  Purpose
C  =======
C
C  DLAR2V applies a vector of real plane rotations from both sides to
C  a sequence of 2-by-2 real symmetric matrices, defined by the elements
C  of the vectors x, y and z. For i = 1,2,...,n
C
C     ( x(i)  z(i) ) := (  c(i)  s(i) ) ( x(i)  z(i) ) ( c(i) -s(i) )
C     ( z(i)  y(i) )    ( -s(i)  c(i) ) ( z(i)  y(i) ) ( s(i)  c(i) )
C
C  Arguments
C  =========
C
C  N       (input) INTEGER
C          The number of plane rotations to be applied.
C
C  X       (input/output) DOUBLE PRECISION array,
C                         dimension (1+(N-1)*INCX)
C          The vector x.
C
C  Y       (input/output) DOUBLE PRECISION array,
C                         dimension (1+(N-1)*INCX)
C          The vector y.
C
C  Z       (input/output) DOUBLE PRECISION array,
C                         dimension (1+(N-1)*INCX)
C          The vector z.
C
C  INCX    (input) INTEGER
C          The increment between elements of X, Y and Z. INCX > 0.
C
C  C       (input) DOUBLE PRECISION array, dimension (1+(N-1)*INCC)
C          The cosines of the plane rotations.
C
C  S       (input) DOUBLE PRECISION array, dimension (1+(N-1)*INCC)
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
      INTEGER           INCC, INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION  C(*), S(*), X(*), Y(*), Z(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  CI, SI, T1, T2, T3, T4, T5, T6, XI, YI, ZI
      INTEGER           I, IC, IX
C     .. Executable Statements ..
C
      IX = 1
      IC = 1
      DO 20 I = 1, N
         XI = X(IX)
         YI = Y(IX)
         ZI = Z(IX)
         CI = C(IC)
         SI = S(IC)
         T1 = SI*ZI
         T2 = CI*ZI
         T3 = T2 - SI*XI
         T4 = T2 + SI*YI
         T5 = CI*XI + T1
         T6 = CI*YI - T1
         X(IX) = CI*T5 + SI*T4
         Y(IX) = CI*T6 - SI*T3
         Z(IX) = CI*T4 - SI*T5
         IX = IX + INCX
         IC = IC + INCC
   20 CONTINUE
C
C     End of F08HEX (DLAR2V)
C
      RETURN
      END
