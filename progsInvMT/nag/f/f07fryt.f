      SUBROUTINE F07FRY(N,X,INCX)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             ZLACGV(N,X,INCX)
C
C  Purpose
C  =======
C
C     ZLACGV conjugates a complex vector of length N.
C
C  Arguments
C  =========
C
C  N      - INTEGER
C           On entry, N specifies the length of the vector X.
C           Unchanged on exit.
C
C  X      - COMPLEX array, dimension( N )
C           On entry, X contains the vector to be conjugated.
C           On return, X contains conjg(X).
C
C  INCX   - INTEGER
C           On entry, INCX specifies the spacing betweeen successive
C           elements of X.
C           Unchanged on exit.
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
C     Courant Institute, NAG Ltd., and Rice University
C
C     .. Scalar Arguments ..
      INTEGER           INCX, N
C     .. Array Arguments ..
      COMPLEX*16        X(*)
C     .. Local Scalars ..
      INTEGER           I, IOFF
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG
C     .. Executable Statements ..
C
      IF (INCX.EQ.1) THEN
         DO 20 I = 1, N
            X(I) = DCONJG(X(I))
   20    CONTINUE
      ELSE
         IOFF = 1
         IF (INCX.LT.0) IOFF = 1 - (N-1)*INCX
         DO 40 I = 1, N
            X(IOFF) = DCONJG(X(IOFF))
            IOFF = IOFF + INCX
   40    CONTINUE
      END IF
      RETURN
C
C     End of F07FRY (ZLACGV)
C
      END
