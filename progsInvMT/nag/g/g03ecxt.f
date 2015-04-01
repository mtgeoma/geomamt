      SUBROUTINE G03ECX(N,D,INC,ITH,JTH,DMIN)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Selects the two clusters to merge
C
C     Only units with IND gt 0 are considered
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DMIN
      INTEGER           ITH, JTH, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(N*(N-1)/2)
      INTEGER           INC(N)
C     .. Local Scalars ..
      INTEGER           I, J, K, L
C     .. Executable Statements ..
C
      DO 20 I = 2, N
         IF (INC(I).GT.0) THEN
            L = I
            GO TO 40
         END IF
   20 CONTINUE
   40 CONTINUE
      K = (L-1)*(L-2)/2 + 1
      DMIN = D(K)
      ITH = L
      JTH = 1
      DO 80 I = L, N
         IF (INC(I).GT.0) THEN
            DO 60 J = 1, I - 1
               IF (INC(J).GT.0) THEN
                  IF (D(K).LE.DMIN) THEN
                     DMIN = D(K)
                     ITH = I
                     JTH = J
                  END IF
               END IF
               K = K + 1
   60       CONTINUE
         ELSE
            K = K + I - 1
         END IF
   80 CONTINUE
      RETURN
      END
