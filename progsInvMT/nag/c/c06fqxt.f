      SUBROUTINE C06FQX(A,B,M,N,Q,NQ,TRIG)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     Hermitian to Real Fast Fourier Transform Kernel Driver
C
C     Self-sorting, decimation in frequency
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
      EXTERNAL          C06FQQ, C06FQR, C06FQS, C06FQT, C06FQU, C06FQV,
     *                  C06FQW
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT, DBLE
C     .. Executable Statements ..
      INA = .TRUE.
      P = 1
      R = N
      IF (N.EQ.1) RETURN
      CALL C06FQQ(A,M,N)
      DO 20 I = 1, NQ
         QI = Q(I)
         R = R/QI
         IF (INA) THEN
            IF (QI.EQ.2) THEN
               CALL C06FQW(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.3) THEN
               CALL C06FQV(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.4) THEN
               CALL C06FQU(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.5) THEN
               CALL C06FQT(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.6) THEN
               CALL C06FQS(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE
               CALL C06FQR(A,B,M*P,QI,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            END IF
         ELSE
            IF (QI.EQ.2) THEN
               CALL C06FQW(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.3) THEN
               CALL C06FQV(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.4) THEN
               CALL C06FQU(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.5) THEN
               CALL C06FQT(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.6) THEN
               CALL C06FQS(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE
               CALL C06FQR(B,A,M*P,QI,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            END IF
         END IF
         INA = .NOT. INA
         P = P*QI
   20 CONTINUE
C
      FACTOR = 2.0D0/SQRT(DBLE(N))
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
