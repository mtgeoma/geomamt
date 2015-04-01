      SUBROUTINE C06FKF(JOB,X,Y,N,WORK,IFAIL)
C     MARK 11 RELEASE. NAG COPYRIGHT 1984.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     IF JOB = 1, CONVOLUTION OF 2 REAL VECTORS.
C     IF JOB = 2, CORRELATION OF 2 REAL VECTORS.
C
C     (USING WORK ARRAY FOR EXTRA SPEED)
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06FKF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, JOB, N
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(N), X(N), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  SQRTN, XI, XN, XR, YI, YR
      INTEGER           I, IERROR, J, ND2, NJ, PMAX, TWOGRP
C     .. Local Arrays ..
      INTEGER           RFACT(21), TFACT(21)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C06EAW, C06EBW, C06FAY, C06FAZ
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Data statements ..
      DATA              PMAX/19/
      DATA              TWOGRP/8/
C     .. Executable Statements ..
      IF (JOB.LT.1 .OR. JOB.GT.2) GO TO 140
      IF (N.LE.1) GO TO 160
C
C     COMPUTE DFT OF X AND Y
C
      CALL C06FAZ(N,PMAX,TWOGRP,TFACT,RFACT,IERROR)
      IF (IERROR.NE.0) GO TO 180
      CALL C06FAY(X,N,RFACT,WORK)
      CALL C06EAW(X,N,TFACT)
      CALL C06FAY(Y,N,RFACT,WORK)
      CALL C06EAW(Y,N,TFACT)
C
C     COMPUTE DFT OF CONVOLUTION OR CORRELATION
C
      ND2 = (N+1)/2
      IF (ND2.LE.1) GO TO 80
      IF (JOB.EQ.2) GO TO 40
      DO 20 J = 2, ND2
         NJ = N + 2 - J
         XR = X(J)
         XI = X(NJ)
         YR = Y(J)
         YI = Y(NJ)
         Y(J) = XR*YR - XI*YI
         Y(NJ) = XR*YI + XI*YR
         X(J) = Y(J)
         X(NJ) = -Y(NJ)
   20 CONTINUE
      GO TO 80
   40 DO 60 J = 2, ND2
         NJ = N + 2 - J
         XR = X(J)
         XI = X(NJ)
         YR = Y(J)
         YI = Y(NJ)
         Y(J) = XR*YR + XI*YI
         Y(NJ) = XR*YI - XI*YR
         X(J) = Y(J)
         X(NJ) = -Y(NJ)
   60 CONTINUE
   80 Y(1) = X(1)*Y(1)
      X(1) = Y(1)
      IF (ND2*2.NE.N) GO TO 100
      Y(ND2+1) = X(ND2+1)*Y(ND2+1)
      X(ND2+1) = Y(ND2+1)
C
C     COMPUTE CONVOLUTION OR CORRELATION BY INVERSE DFT
C
  100 CALL C06EBW(X,N,TFACT)
      CALL C06FAY(X,N,RFACT,WORK)
C
C     SCALE
C
      XN = N
      SQRTN = SQRT(XN)
      DO 120 I = 1, N
         X(I) = X(I)/XN
         Y(I) = Y(I)/SQRTN
  120 CONTINUE
C
      IFAIL = 0
      RETURN
C
  140 IERROR = 4
      GO TO 180
  160 IERROR = 3
  180 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
