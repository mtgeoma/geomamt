      SUBROUTINE C06ECW(X,Y,PTS,FACTOR)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     COMPLEX FOURIER TRANSFORM KERNEL DRIVER
C     .. Scalar Arguments ..
      INTEGER           PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X(PTS), Y(PTS)
      INTEGER           FACTOR(21)
C     .. Local Scalars ..
      INTEGER           F, M, M1, M2, M3, M4, M5, M6, M7, P
C     .. External Subroutines ..
      EXTERNAL          C06ECQ, C06ECR, C06ECS, C06ECT, C06ECU, C06ECV
C     .. Executable Statements ..
      F = 0
      M = PTS
   20 CONTINUE
      F = F + 1
      P = FACTOR(F)
      IF (P.EQ.0) RETURN
      IF (P.EQ.1) GO TO 20
      M = M/P
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
      CALL C06ECV(X(1),Y(1),PTS,X(M+1),Y(M+1),M1,M)
      GO TO 20
C
   60 CONTINUE
      CALL C06ECU(X(1),Y(1),PTS,X(M+1),Y(M+1),M1,X(2*M+1),Y(2*M+1),M2,M)
      GO TO 20
C
   80 CONTINUE
      CALL C06ECT(X(1),Y(1),PTS,X(M+1),Y(M+1),M1,X(2*M+1),Y(2*M+1)
     *            ,M2,X(3*M+1),Y(3*M+1),M3,M)
      GO TO 20
C
  100 CONTINUE
      CALL C06ECS(X(1),Y(1),PTS,X(M+1),Y(M+1),M1,X(2*M+1),Y(2*M+1)
     *            ,M2,X(3*M+1),Y(3*M+1),M3,X(4*M+1),Y(4*M+1),M4,M)
      GO TO 20
C
  120 CONTINUE
      CALL C06ECR(X(1),Y(1),PTS,X(M+1),Y(M+1),M1,X(2*M+1),Y(2*M+1)
     *            ,M2,X(3*M+1),Y(3*M+1),M3,X(4*M+1),Y(4*M+1),M4,X(5*M+1)
     *            ,Y(5*M+1),M5,X(6*M+1),Y(6*M+1),M6,X(7*M+1),Y(7*M+1)
     *            ,M7,M)
      GO TO 20
C
  140 CONTINUE
      CALL C06ECQ(X,Y,PTS,M,P)
      GO TO 20
C
      END
