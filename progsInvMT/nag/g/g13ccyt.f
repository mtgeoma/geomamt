      SUBROUTINE G13CCY(X,K,M,S1,S2)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G13CCY RETRIEVES THE REAL AND IMAGINARY PARTS
C     OF Z  FROM THE ARRAY X.
C
C     X  - ARRAY HOLDING HERMITIAN SEQUENCE OF COMPLEX NOS.
C     K  - LENGTH OF ARRAY X
C     M  - SUBSCRIPT OF Z
C     S1 - TO RETURN REAL PART
C     S2 - TO RETURN IMAGINARY PART
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  S1, S2
      INTEGER           K, M
C     .. Array Arguments ..
      DOUBLE PRECISION  X(K)
C     .. Local Scalars ..
      INTEGER           M1, M2, M3, N
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Executable Statements ..
      M1 = ABS(M) + 1
      N = K + 2 - M1
      M2 = MIN(M1,N)
      M3 = MAX(M1,N)
      S1 = X(M2)
      S2 = 0.0D0
      IF (M3.NE.(K+1) .AND. M1.NE.N) S2 = X(M3)
      IF (M.LT.0) S2 = -S2
      IF (M1.GT.N) S2 = -S2
      RETURN
      END
