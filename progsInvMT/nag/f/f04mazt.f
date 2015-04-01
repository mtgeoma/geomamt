      SUBROUTINE F04MAZ(N,NZ,A,INI,INJ,IAF,AF,DF,INJF,IK,B,R,E,F,G,KMAX,
     *                  EPS,LROW,IERR)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS
      INTEGER           IAF, IERR, KMAX, LROW, N, NZ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NZ), AF(IAF), B(N), DF(N), E(N), F(N), G(N),
     *                  R(N)
      INTEGER           IK(N,2), INI(NZ), INJ(NZ), INJF(IAF)
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION  BB, D0, D1, L, R1, ZERO
      INTEGER           I, KITR
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, DDOT
      EXTERNAL          DNRM2, DDOT
C     .. External Subroutines ..
      EXTERNAL          F06FGF, DSCAL, DCOPY, DAXPY, F04MAX, F04MAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AX02ZA/WMACH
C     .. Save statements ..
      SAVE              /AX02ZA/
C     .. Data statements ..
      DATA              ZERO/0.0D+0/
C     .. Executable Statements ..
C
C     THIS SUBROUTINE PERFORMS THE ITERATIVE PROCEDURE.
C     THE PRECONDITIONED CONJUGATE GRADIENT METHOD IS USED.
C
C
C     COMPUTE THE INITIAL SOLUTION.
C
      IERR = 0
      CALL DCOPY(N,B,1,E,1)
C
      CALL F04MAY(N,AF,INJF,IAF,DF,IK,E,LROW)
C
C     COMPUTE THE RESIDUALS AND INSERT THE INITIAL SOLUTION IN B.
C
      CALL F04MAX(A,INI,INJ,NZ,N,E,R)
C
      DO 20 I = 1, N
         R(I) = R(I) - B(I)
   20 CONTINUE
C
      R1 = DNRM2(N,R,1)
      CALL DCOPY(N,R,1,G,1)
      CALL DCOPY(N,E,1,B,1)
C
      KITR = 0
      IF (R1.LE.EPS) GO TO 100
C
C     INITIALIZE E AND G.
C
      CALL F04MAY(N,AF,INJF,IAF,DF,IK,G,LROW)
C
      CALL DCOPY(N,G,1,E,1)
      CALL F06FGF(N,E,1)
      D0 = DDOT(N,R,1,G,1)
C
C     START ITERATION LOOP
C
   40 KITR = KITR + 1
C
      CALL F04MAX(A,INI,INJ,NZ,N,E,F)
C
      L = DDOT(N,E,1,F,1)
      IF (ABS(L).GT.ABS(D0)*WMACH(3)) GO TO 60
      KMAX = KITR
      IERR = 3
      RETURN
   60 L = D0/L
C
C     ADJUST B, G AND R.
C
      CALL DAXPY(N,L,E,1,B,1)
      CALL DAXPY(N,L,F,1,R,1)
      R1 = DNRM2(N,R,1)
      CALL DCOPY(N,R,1,G,1)
C
C     CONTROL THE RESIDUAL.
C
      IF (R1.LE.EPS .OR. KITR.GE.KMAX) GO TO 100
C
C     PROCEED ITERATION .
C
      CALL F04MAY(N,AF,INJF,IAF,DF,IK,G,LROW)
C
      D1 = DDOT(N,R,1,G,1)
C
      BB = D1/D0
      D0 = D1
C
      CALL DSCAL(N,BB,E,1)
C
      DO 80 I = 1, N
         E(I) = E(I) - G(I)
   80 CONTINUE
      GO TO 40
C
C     ITERATION LOOP TERMINATES.
C
  100 KMAX = KITR
      EPS = R1
      RETURN
C
C     END OF F04MAZ. ( MA31F. )
C
      END
