      SUBROUTINE F08HSX(N,X,Y,Z,INCX,C,S,INCC)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZLAR2V(N,X,Y,Z,INCX,C,S,INCC)
C
C  Purpose
C  =======
C
C  ZLAR2V applies a vector of complex plane rotations with real cosines
C  from both sides to a sequence of 2-by-2 complex Hermitian matrices,
C  defined by the elements of the vectors x, y and z. For i = 1,2,...,n
C
C     (       x(i)  z(i) ) :=
C     ( conjg(z(i)) y(i) )
C
C       (  c(i) conjg(s(i)) ) (       x(i)  z(i) ) ( c(i) -conjg(s(i)) )
C       ( -s(i)       c(i)  ) ( conjg(z(i)) y(i) ) ( s(i)        c(i)  )
C
C  Arguments
C  =========
C
C  N       (input) INTEGER
C          The number of plane rotations to be applied.
C
C  X       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCX)
C          The vector x; the elements of x are assumed to be real.
C
C  Y       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCX)
C          The vector y; the elements of y are assumed to be real.
C
C  Z       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCX)
C          The vector z.
C
C  INCX    (input) INTEGER
C          The increment between elements of X, Y and Z. INCX > 0.
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
      INTEGER           INCC, INCX, N
C     .. Array Arguments ..
      COMPLEX*16        S(*), X(*), Y(*), Z(*)
      DOUBLE PRECISION  C(*)
C     .. Local Scalars ..
      COMPLEX*16        SI, T2, T3, T4, ZI
      DOUBLE PRECISION  CI, SII, SIR, T1I, T1R, T5, T6, XI, YI, ZII, ZIR
      INTEGER           I, IC, IX
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, DCMPLX, DCONJG, DIMAG
C     .. Executable Statements ..
C
      IX = 1
      IC = 1
      DO 20 I = 1, N
         XI = DBLE(X(IX))
         YI = DBLE(Y(IX))
         ZI = Z(IX)
         ZIR = DBLE(ZI)
         ZII = DIMAG(ZI)
         CI = C(IC)
         SI = S(IC)
         SIR = DBLE(SI)
         SII = DIMAG(SI)
         T1R = SIR*ZIR - SII*ZII
         T1I = SIR*ZII + SII*ZIR
         T2 = CI*ZI
         T3 = T2 - DCONJG(SI)*XI
         T4 = DCONJG(T2) + SI*YI
         T5 = CI*XI + T1R
         T6 = CI*YI - T1R
         X(IX) = CI*T5 + (SIR*DBLE(T4)+SII*DIMAG(T4))
         Y(IX) = CI*T6 - (SIR*DBLE(T3)-SII*DIMAG(T3))
         Z(IX) = CI*T3 + DCONJG(SI)*DCMPLX(T6,T1I)
         IX = IX + INCX
         IC = IC + INCC
   20 CONTINUE
      RETURN
C
C     End of F08HSX (ZLAR2V)
C
      END
