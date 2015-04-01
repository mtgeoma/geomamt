      SUBROUTINE C06EAW(X,PTS,FACTOR)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     REAL FOURIER TRANSFORM KERNEL DRIVER
C     .. Scalar Arguments ..
      INTEGER           PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X(PTS)
      INTEGER           FACTOR(21)
C     .. Local Scalars ..
      INTEGER           F, M, M1, M2, M3, M4, M5, M6, M7, P
C     .. External Subroutines ..
      EXTERNAL          C06EAQ, C06EAR, C06EAS, C06EAT, C06EAU, C06EAV
C     .. Executable Statements ..
      M = 1
      P = 1
      F = 0
   20 CONTINUE
      M = M*P
      F = F + 1
      P = FACTOR(F)
      IF (P.EQ.0) RETURN
      IF (P.EQ.1) GO TO 20
      M1 = PTS - M
      M2 = M1 - M
      M3 = M2 - M
      M4 = M3 - M
      M5 = M4 - M
      M6 = M5 - M
      M7 = M6 - M
      IF (P.EQ.2) GO TO 40
      IF (P.EQ.3) GO TO 60
      IF (P.EQ.4) GO TO 80
      IF (P.EQ.5) GO TO 100
      IF (P.EQ.8) GO TO 120
      GO TO 140
C
   40 CONTINUE
      CALL C06EAV(X(1),PTS,X(M+1),M1,M)
      GO TO 20
C
   60 CONTINUE
      CALL C06EAU(X(1),PTS,X(M+1),M1,X(2*M+1),M2,M)
      GO TO 20
C
   80 CONTINUE
      CALL C06EAT(X(1),PTS,X(M+1),M1,X(2*M+1),M2,X(3*M+1),M3,M)
      GO TO 20
C
  100 CONTINUE
      CALL C06EAS(X(1),PTS,X(M+1),M1,X(2*M+1),M2,X(3*M+1),M3,X(4*M+1)
     *            ,M4,M)
      GO TO 20
C
  120 CONTINUE
      CALL C06EAR(X(1),PTS,X(M+1),M1,X(2*M+1),M2,X(3*M+1),M3,X(4*M+1)
     *            ,M4,X(5*M+1),M5,X(6*M+1),M6,X(7*M+1),M7,M)
      GO TO 20
C
  140 CONTINUE
      CALL C06EAQ(X,PTS,M,P)
      GO TO 20
C
      END
