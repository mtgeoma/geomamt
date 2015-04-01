      SUBROUTINE C06FPX(A,B,M,N,Q,NQ,TRIG)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     Real to Hermitian Fast Fourier Transform Kernel Driver
C
C     Mixed-radix, self-sorting, decimation in time
C
C     .. Scalar Arguments ..
      INTEGER           M, N, NQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:M*N-1), B(0:M*N-1), TRIG(0:2*N-1)
      INTEGER           Q(NQ)
C     .. Local Scalars ..
      DOUBLE PRECISION  FACTOR
      INTEGER           I, P, QI, R
      LOGICAL           INA
C     .. External Subroutines ..
      EXTERNAL          C06FPR, C06FPS, C06FPT, C06FPU, C06FPV, C06FPW
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT, DBLE
C     .. Executable Statements ..
      INA = .TRUE.
      P = N
      R = 1
      IF (N.EQ.1) RETURN
      DO 20 I = NQ, 1, -1
         QI = Q(I)
         P = P/QI
         IF (INA) THEN
            IF (QI.EQ.2) THEN
               CALL C06FPW(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.3) THEN
               CALL C06FPV(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.4) THEN
               CALL C06FPU(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.5) THEN
               CALL C06FPT(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.6) THEN
               CALL C06FPS(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE
               CALL C06FPR(A,B,M*P,QI,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            END IF
         ELSE
            IF (QI.EQ.2) THEN
               CALL C06FPW(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.3) THEN
               CALL C06FPV(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.4) THEN
               CALL C06FPU(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.5) THEN
               CALL C06FPT(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.6) THEN
               CALL C06FPS(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE
               CALL C06FPR(B,A,M*P,QI,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            END IF
         END IF
         INA = .NOT. INA
         R = R*QI
   20 CONTINUE
C
      FACTOR = 1.0D0/SQRT(DBLE(N))
      IF (INA) THEN
         DO 40 I = 0, M*N - 1
            A(I) = A(I)*FACTOR
   40    CONTINUE
      ELSE
         DO 60 I = 0, M*N - 1
            A(I) = B(I)*FACTOR
   60    CONTINUE
      END IF
C
      RETURN
      END
