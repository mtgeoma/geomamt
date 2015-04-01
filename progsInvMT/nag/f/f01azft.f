      SUBROUTINE F01AZF(N,M1,M2,A,IA,Z,IZI)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C     MARK 5C REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     TRBAK3
C     THIS SUBROUTINE PERFORMS, ON THE MATRIX OF EIGENVECTORS, Z,
C     STORED
C     IN COLUMNS M1 TO M2 OF THE ARRAY Z(N,N), A BACKTRANSFORMATION
C     TO
C     FORM THE EIGENVECTORS OF THE ORIGINAL SYMMETRIC MATRIX FROM
C     THE EIGENVECTORS OF THE TRIDIAGONAL MATRIX. THE NEW VECTORS
C     ARE OVERWRITTEN ON THE OLD ONES. THE DETAILS OF THE
C     HOUSEHOLDER REDUCTION MUST BE STORED IN THE ARRAY
C     A(N(N+1)/2), AS LEFT BY THE SUBROUTINE TRED3, F01AYF. IF Z
C     DENOTES ANY COLUMN OF THE RESULTANT MATRIX Z, THEN Z
C     SATISFIES ZT*Z = Z(INPUT)T * Z(INPUT).
C     IST. MARCH  1972
C
C     .. Scalar Arguments ..
      INTEGER           IA, IZI, M1, M2, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA), Z(IZI,M2)
C     .. Local Scalars ..
      DOUBLE PRECISION  H, S
      INTEGER           I, IPOS, IZ, J, K, L
C     .. Executable Statements ..
      IF (N.EQ.1) RETURN
      DO 80 I = 2, N
         L = I - 1
         IZ = I*L/2
         IPOS = IZ + I
         H = A(IPOS)
         IF (H.EQ.0.0D0) GO TO 80
         DO 60 J = M1, M2
            S = 0.0D0
            DO 20 K = 1, L
               IPOS = IZ + K
               S = S + A(IPOS)*Z(K,J)
   20       CONTINUE
            S = S/H
            DO 40 K = 1, L
               IPOS = IZ + K
               Z(K,J) = Z(K,J) - S*A(IPOS)
   40       CONTINUE
   60    CONTINUE
   80 CONTINUE
      RETURN
      END