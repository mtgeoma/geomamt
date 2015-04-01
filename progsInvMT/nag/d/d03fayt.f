      SUBROUTINE D03FAY(LBDCND,L,MBDCND,M,N,A,B,C,LDIMF,MDIMF,F,XRT,YRT,
     *                  C1,C2,WORK1,JWORK1,WORK2,JWORK2,TRIGL,JTRIGL,
     *                  TRIGM,JTRIGM)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  C1, C2
      INTEGER           JTRIGL, JTRIGM, JWORK1, JWORK2, L, LBDCND,
     *                  LDIMF, M, MBDCND, MDIMF, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N), B(N), C(N), F(LDIMF,MDIMF,*),
     *                  TRIGL(JTRIGL), TRIGM(JTRIGM), WORK1(JWORK1),
     *                  WORK2(JWORK2), XRT(L), YRT(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  DI, DJ, DX, DY, PI
      INTEGER           I, IFAIL, J, K
      CHARACTER*1       INIT
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. External Subroutines ..
      EXTERNAL          C06FPF, C06FQF, C06GQF, C06HAF, C06HBF, C06HCF,
     *                  C06HDF, D03FAX
C     .. Intrinsic Functions ..
      INTRINSIC         MOD, DBLE, SIN
C     .. Executable Statements ..
      PI = X01AAF(PI)
C
C     Generate scale factors
C
      IF (LBDCND.EQ.0) THEN
         DX = 0.5D0*PI/DBLE(L)
      ELSE IF (LBDCND.EQ.1) THEN
         DX = 0.5D0*PI/DBLE(L+1)
      ELSE IF (LBDCND.EQ.2) THEN
         DX = 0.5D0*PI/DBLE(L)
      ELSE IF (LBDCND.EQ.3) THEN
         DX = 0.5D0*PI/DBLE(L-1)
      ELSE IF (LBDCND.EQ.4) THEN
         DX = 0.5D0*PI/DBLE(L)
      END IF
      IF (LBDCND.EQ.0) THEN
C
C        Here are the equations for the periodic
C        boundary conditions in x
C
         XRT(1) = 0.0D0
         DO 20 I = 2, (L+1)/2
            XRT(I) = -4.0D0*C1*(SIN(DBLE(2*(I-1))*DX))**2
            XRT(L+2-I) = XRT(I)
   20    CONTINUE
         IF (MOD(L,2).EQ.0) XRT((L/2)+1) = -4.0D0*C1
      ELSE
         IF (LBDCND.EQ.1) DI = 0.0D0
         IF (LBDCND.EQ.2 .OR. LBDCND.EQ.4) DI = 0.5D0
         IF (LBDCND.EQ.3) DI = 1.0D0
         DO 40 I = 1, L
            XRT(I) = -4.0D0*C1*(SIN((DBLE(I)-DI)*DX))**2
   40    CONTINUE
      END IF
      IF (MBDCND.EQ.0) THEN
         DY = 0.5D0*PI/DBLE(M)
      ELSE IF (MBDCND.EQ.1) THEN
         DY = 0.5D0*PI/DBLE(M+1)
      ELSE IF (MBDCND.EQ.2) THEN
         DY = 0.5D0*PI/DBLE(M)
      ELSE IF (MBDCND.EQ.3) THEN
         DY = 0.5D0*PI/DBLE(M-1)
      ELSE IF (MBDCND.EQ.4) THEN
         DY = 0.5D0*PI/DBLE(M)
      END IF
      IF (MBDCND.EQ.0) THEN
C
C        Here are the equations for the periodic
C        boundary conditions in Y
C
         YRT(1) = 0.0D0
         DO 60 I = 2, (M+1)/2
            YRT(I) = -4.0D0*C2*(SIN(DBLE(2*(I-1))*DY))**2
            YRT(M+2-I) = YRT(I)
   60    CONTINUE
         IF (MOD(M,2).EQ.0) YRT((M/2)+1) = -4.0D0*C2
      ELSE
         IF (MBDCND.EQ.1) DJ = 0.0D0
         IF (MBDCND.EQ.2 .OR. MBDCND.EQ.4) DJ = 0.5D0
         IF (MBDCND.EQ.3) DJ = 1.0D0
         DO 80 J = 1, M
            YRT(J) = -4.0D0*C2*(SIN((DBLE(J)-DJ)*DY))**2
   80    CONTINUE
      END IF
C
C     Transform X
C
      DO 180 J = 1, M
         INIT = 'Subsequent'
         IF (J.EQ.1) INIT = 'Initial'
         DO 120 K = 1, N
            DO 100 I = 1, L
               WORK1((I-1)*N+K) = F(I,J,K)
  100       CONTINUE
  120    CONTINUE
         IFAIL = 0
         IF (LBDCND.EQ.0) THEN
            CALL C06FPF(N,L,WORK1,INIT,TRIGL,WORK2,IFAIL)
         ELSE IF (LBDCND.EQ.1) THEN
            CALL C06HAF(N,L+1,WORK1,INIT,TRIGL,WORK2,IFAIL)
         ELSE IF (LBDCND.EQ.2) THEN
            CALL C06HCF('Forward',N,L,WORK1,INIT,TRIGL,WORK2,IFAIL)
         ELSE IF (LBDCND.EQ.3) THEN
            CALL C06HBF(N,L-1,WORK1,INIT,TRIGL,WORK2,IFAIL)
         ELSE IF (LBDCND.EQ.4) THEN
            CALL C06HDF('Forward',N,L,WORK1,INIT,TRIGL,WORK2,IFAIL)
         END IF
         DO 160 K = 1, N
            DO 140 I = 1, L
               F(I,J,K) = WORK1((I-1)*N+K)
  140       CONTINUE
  160    CONTINUE
  180 CONTINUE
C
C     Transform Y
C
      DO 280 I = 1, L
         INIT = 'Subsequent'
         IF (I.EQ.1) INIT = 'Initial'
         DO 220 K = 1, N
            DO 200 J = 1, M
               WORK1((J-1)*N+K) = F(I,J,K)
  200       CONTINUE
  220    CONTINUE
         IFAIL = 0
         IF (MBDCND.EQ.0) THEN
            CALL C06FPF(N,M,WORK1,INIT,TRIGM,WORK2,IFAIL)
         ELSE IF (MBDCND.EQ.1) THEN
            CALL C06HAF(N,M+1,WORK1,INIT,TRIGM,WORK2,IFAIL)
         ELSE IF (MBDCND.EQ.2) THEN
            CALL C06HCF('Forward',N,M,WORK1,INIT,TRIGM,WORK2,IFAIL)
         ELSE IF (MBDCND.EQ.3) THEN
            CALL C06HBF(N,M-1,WORK1,INIT,TRIGM,WORK2,IFAIL)
         ELSE IF (MBDCND.EQ.4) THEN
            CALL C06HDF('Forward',N,M,WORK1,INIT,TRIGM,WORK2,IFAIL)
         END IF
         DO 260 K = 1, N
            DO 240 J = 1, M
               F(I,J,K) = WORK1((J-1)*N+K)
  240       CONTINUE
  260    CONTINUE
  280 CONTINUE
C
C     Solve tridiagonal systems in Z
C
      DO 380 I = 1, L
C
         DO 320 K = 1, N
            DO 300 J = 1, M
               WORK1(J+(K-1)*M) = F(I,J,K)
  300       CONTINUE
  320    CONTINUE
         CALL D03FAX(M,N,A,B,C,WORK1,WORK2,XRT(I),YRT)
         DO 360 K = 1, N
            DO 340 J = 1, M
               F(I,J,K) = WORK1(J+(K-1)*M)
  340       CONTINUE
  360    CONTINUE
C
  380 CONTINUE
C
C     Un-transform Y
C
      DO 480 I = 1, L
         DO 420 K = 1, N
            DO 400 J = 1, M
               WORK1((J-1)*N+K) = F(I,J,K)
  400       CONTINUE
  420    CONTINUE
         IFAIL = 0
         IF (MBDCND.EQ.0) THEN
            CALL C06GQF(N,M,WORK1,IFAIL)
            CALL C06FQF(N,M,WORK1,'Subsequent',TRIGM,WORK2,IFAIL)
         ELSE IF (MBDCND.EQ.1) THEN
            CALL C06HAF(N,M+1,WORK1,'Subsequent',TRIGM,WORK2,IFAIL)
         ELSE IF (MBDCND.EQ.2) THEN
            CALL C06HCF('Backwards',N,M,WORK1,'Subsequent',TRIGM,WORK2,
     *                  IFAIL)
         ELSE IF (MBDCND.EQ.3) THEN
            CALL C06HBF(N,M-1,WORK1,'Subsequent',TRIGM,WORK2,IFAIL)
         ELSE IF (MBDCND.EQ.4) THEN
            CALL C06HDF('Backward',N,M,WORK1,'Subsequent',TRIGM,WORK2,
     *                  IFAIL)
         END IF
         DO 460 K = 1, N
            DO 440 J = 1, M
               F(I,J,K) = WORK1((J-1)*N+K)
  440       CONTINUE
  460    CONTINUE
  480 CONTINUE
C
C     Un-transform X
C
      DO 580 J = 1, M
         DO 520 K = 1, N
            DO 500 I = 1, L
               WORK1((I-1)*N+K) = F(I,J,K)
  500       CONTINUE
  520    CONTINUE
         IFAIL = 0
         IF (LBDCND.EQ.0) THEN
            CALL C06GQF(N,L,WORK1,IFAIL)
            CALL C06FQF(N,L,WORK1,'Subsequent',TRIGL,WORK2,IFAIL)
         ELSE IF (LBDCND.EQ.1) THEN
            CALL C06HAF(N,L+1,WORK1,'Subsequent',TRIGL,WORK2,IFAIL)
         ELSE IF (LBDCND.EQ.2) THEN
            CALL C06HCF('Backward',N,L,WORK1,'Subsequent',TRIGL,WORK2,
     *                  IFAIL)
         ELSE IF (LBDCND.EQ.3) THEN
            CALL C06HBF(N,L-1,WORK1,'Subsequent',TRIGL,WORK2,IFAIL)
         ELSE IF (LBDCND.EQ.4) THEN
            CALL C06HDF('Backward',N,L,WORK1,'Subsequent',TRIGL,WORK2,
     *                  IFAIL)
         END IF
         DO 560 K = 1, N
            DO 540 I = 1, L
               F(I,J,K) = WORK1((I-1)*N+K)
  540       CONTINUE
  560    CONTINUE
  580 CONTINUE
C
      RETURN
      END
