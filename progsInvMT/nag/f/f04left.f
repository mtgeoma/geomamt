      SUBROUTINE F04LEF(JOB,N,A,B,C,D,IN,Y,TOL,IFAIL)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-743 (DEC 1989).
C     MARK 16A REVISED. IER-1006 (JUN 1993).
C
C     F04LEF SOLVES ONE OF THE SYSTEMS OF EQUATIONS
C
C     T*X = Y ,   ( T**T )*X = Y ,   U*X = Y ,
C
C     WHERE T IS AN N BY N TRIDIAGONAL MATRIX THAT HAS BEEN
C     FACTORIZED AS
C
C     T = P*L*U ,
C
C     BY ROUTINE F01LEF.
C
C     FOR A DESCRIPTION OF THE PARAMETERS AND USE OF THIS ROUTINE SEE
C     THE NAG LIBRARY MANUAL.
C
C     -- WRITTEN ON 19-JANUARY-1983.  S.J.HAMMARLING.
C
C     NAG FORTRAN 66 GENERAL PURPOSE ROUTINE.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04LEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL
      INTEGER           IFAIL, JOB, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*), B(*), C(*), D(*), Y(*)
      INTEGER           IN(*)
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION  AK, PERT, TEMP, TWO, ZERO
      INTEGER           IERR, K, KK
      LOGICAL           FAIL
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  F06BLF
      INTEGER           P01ABF
      EXTERNAL          F06BLF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          X02ZAZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Common blocks ..
      COMMON            /AX02ZA/WMACH
C     .. Save statements ..
      SAVE              /AX02ZA/
C     .. Data statements ..
      DATA              TWO/2.0D+0/, ZERO/0.0D+0/
C     .. Executable Statements ..
C
      IF (N.GT.0 .AND. ABS(JOB).LT.4 .AND. JOB.NE.0) GO TO 20
      IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
   20 CONTINUE
C
      CALL X02ZAZ
C
      IF (JOB.GT.0) GO TO 100
      IF (TOL.GT.ZERO) GO TO 80
      TOL = ABS(A(1))
      IF (N.GT.1) TOL = MAX(TOL,ABS(A(2)),ABS(B(2)))
      IF (N.LT.3) GO TO 60
      DO 40 K = 3, N
         TOL = MAX(TOL,ABS(A(K)),ABS(B(K)),ABS(D(K)))
   40 CONTINUE
   60 CONTINUE
      TOL = TOL*WMACH(3)
      IF (TOL.EQ.ZERO) TOL = WMACH(3)
   80 CONTINUE
  100 CONTINUE
C
      IF (ABS(JOB).NE.1 .OR. N.EQ.1) GO TO 200
      DO 180 K = 2, N
         IF (IN(K-1).EQ.1) GO TO 120
         Y(K) = Y(K) - C(K)*Y(K-1)
         GO TO 140
  120    CONTINUE
         TEMP = Y(K-1)
         Y(K-1) = Y(K)
         Y(K) = TEMP - C(K)*Y(K)
  140    CONTINUE
  180 CONTINUE
  200 CONTINUE
      IF (JOB.NE.1 .AND. JOB.NE.3) GO TO 280
      Y(N) = F06BLF(Y(N),A(N),FAIL)
      IF (FAIL) GO TO 840
      IF (N.EQ.1) GO TO 260
      Y(N-1) = F06BLF(Y(N-1)-B(N)*Y(N),A(N-1),FAIL)
      IF (FAIL) GO TO 860
      IF (N.EQ.2) GO TO 240
      K = N - 2
      DO 220 KK = 3, N
         Y(K) = F06BLF(Y(K)-B(K+1)*Y(K+1)-D(K+2)*Y(K+2),A(K),FAIL)
         IF (FAIL) GO TO 880
         K = K - 1
  220 CONTINUE
  240 CONTINUE
  260 CONTINUE
      IFAIL = 0
      RETURN
  280 CONTINUE
C
      IF (JOB.NE.(-1) .AND. JOB.NE.(-3)) GO TO 460
      K = N
      DO 440 KK = 1, N
         AK = A(K)
         IF (AK.LT.ZERO) GO TO 300
         PERT = TOL
         GO TO 320
  300    CONTINUE
         PERT = -TOL
  320    CONTINUE
         IF (K.LT.N) GO TO 340
         TEMP = Y(N)
         GO TO 380
  340    IF (K.LT.(N-1)) GO TO 360
         TEMP = Y(N-1) - B(N)*Y(N)
         GO TO 380
  360    CONTINUE
         TEMP = Y(K) - B(K+1)*Y(K+1) - D(K+2)*Y(K+2)
  380    CONTINUE
         Y(K) = F06BLF(TEMP,AK,FAIL)
C        +          WHILE( FAIL ) LOOP
  400    IF ( .NOT. FAIL) GO TO 420
         AK = AK + PERT
         PERT = TWO*PERT
         Y(K) = F06BLF(TEMP,AK,FAIL)
         GO TO 400
C        +          END LOOP
  420    CONTINUE
         K = K - 1
  440 CONTINUE
      IFAIL = 0
      RETURN
  460 CONTINUE
C
      IF (JOB.EQ.(-2)) GO TO 540
      Y(1) = F06BLF(Y(1),A(1),FAIL)
      IF (FAIL) GO TO 900
      IF (N.EQ.1) GO TO 520
      Y(2) = F06BLF(Y(2)-B(2)*Y(1),A(2),FAIL)
      IF (FAIL) GO TO 920
      IF (N.EQ.2) GO TO 500
      DO 480 K = 3, N
         Y(K) = F06BLF(Y(K)-B(K)*Y(K-1)-D(K)*Y(K-2),A(K),FAIL)
         IF (FAIL) GO TO 880
  480 CONTINUE
  500 CONTINUE
  520 CONTINUE
      GO TO 720
  540 CONTINUE
      DO 700 K = 1, N
         AK = A(K)
         IF (AK.LT.ZERO) GO TO 560
         PERT = TOL
         GO TO 580
  560    CONTINUE
         PERT = -TOL
  580    CONTINUE
         IF (K.GT.1) GO TO 600
         TEMP = Y(1)
         GO TO 640
  600    IF (K.GT.2) GO TO 620
         TEMP = Y(2) - B(2)*Y(1)
         GO TO 640
  620    CONTINUE
         TEMP = Y(K) - B(K)*Y(K-1) - D(K)*Y(K-2)
  640    CONTINUE
         Y(K) = F06BLF(TEMP,AK,FAIL)
C        +          WHILE( FAIL ) LOOP
  660    IF ( .NOT. FAIL) GO TO 680
         AK = AK + PERT
         PERT = TWO*PERT
         Y(K) = F06BLF(TEMP,AK,FAIL)
         GO TO 660
C        +          END LOOP
  680    CONTINUE
  700 CONTINUE
  720 CONTINUE
C
      IF (N.EQ.1) GO TO 820
      K = N
      DO 800 KK = 2, N
         IF (IN(K-1).EQ.1) GO TO 740
         Y(K-1) = Y(K-1) - C(K)*Y(K)
         GO TO 760
  740    CONTINUE
         TEMP = Y(K-1)
         Y(K-1) = Y(K)
         Y(K) = TEMP - C(K)*Y(K)
  760    CONTINUE
         K = K - 1
  800 CONTINUE
  820 CONTINUE
      IFAIL = 0
      RETURN
C
  840 IERR = N + 1
      GO TO 940
  860 IERR = N
      GO TO 940
  880 IERR = K + 1
      GO TO 940
  900 IERR = 2
      GO TO 940
  920 IERR = 3
C
  940 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
C
C     END OF F04LEF.
C
      END
