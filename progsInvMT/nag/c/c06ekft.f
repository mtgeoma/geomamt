      SUBROUTINE C06EKF(JOB,X,Y,N,IFAIL)
C     MARK 11 RELEASE. NAG COPYRIGHT 1984.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     IF JOB = 1, CONVOLUTION OF 2 REAL VECTORS.
C     IF JOB = 2, CORRELATION OF 2 REAL VECTORS.
C
C     (NO WORK ARRAY REQUIRED)
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06EKF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, JOB, N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  SQRTN, XI, XN, XR, YI, YR
      INTEGER           I, IERROR, J, ND2, NJ, NSYM, PMAX, TWOGRP
C     .. Local Arrays ..
      INTEGER           FACTOR(21), SYM(21), UNSYM(21)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C06EAW, C06EAY, C06EAZ, C06EBW
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
      CALL C06EAZ(N,PMAX,TWOGRP,FACTOR,SYM,NSYM,UNSYM,IERROR)
      IF (IERROR.NE.0) GO TO 180
      CALL C06EAY(X,N,SYM,NSYM,UNSYM)
      CALL C06EAW(X,N,FACTOR)
      CALL C06EAY(Y,N,SYM,NSYM,UNSYM)
      CALL C06EAW(Y,N,FACTOR)
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
  100 CALL C06EBW(X,N,FACTOR)
      CALL C06EAY(X,N,SYM,NSYM,UNSYM)
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
