      SUBROUTINE C06FRX(X,Y,BR,BI,M,N,Q,NQ,TRIG)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     Multiple complex Fourier transform kernel driver
C
C     Mixed-radix, self-sorting, decimation in frequency
C
C     .. Scalar Arguments ..
      INTEGER           M, N, NQ
C     .. Array Arguments ..
      DOUBLE PRECISION  BI(0:M*N-1), BR(0:M*N-1), TRIG(0:2*N-1),
     *                  X(0:M*N-1), Y(0:M*N-1)
      INTEGER           Q(NQ)
C     .. Local Scalars ..
      DOUBLE PRECISION  FACTOR
      INTEGER           I, P, QI, R
      LOGICAL           INA
C     .. External Subroutines ..
      EXTERNAL          C06FRR, C06FRS, C06FRT, C06FRU, C06FRV, C06FRW
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT, DBLE
C     .. Executable Statements ..
      INA = .TRUE.
      P = 1
      R = N
C
      IF (N.EQ.1) RETURN
      DO 20 I = 1, NQ
         QI = Q(I)
         R = R/QI
         IF (INA) THEN
            IF (QI.EQ.2) THEN
               CALL C06FRW(X,Y,BR,BI,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            ELSE IF (QI.EQ.3) THEN
               CALL C06FRV(X,Y,BR,BI,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            ELSE IF (QI.EQ.4) THEN
               CALL C06FRU(X,Y,BR,BI,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            ELSE IF (QI.EQ.5) THEN
               CALL C06FRT(X,Y,BR,BI,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            ELSE IF (QI.EQ.6) THEN
               CALL C06FRS(X,Y,BR,BI,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            ELSE
               CALL C06FRR(X,Y,BR,BI,M*P,QI,R,TRIG((P-1)*QI*R),
     *                     TRIG(N+(P-1)*QI*R))
            END IF
         ELSE
            IF (QI.EQ.2) THEN
               CALL C06FRW(BR,BI,X,Y,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            ELSE IF (QI.EQ.3) THEN
               CALL C06FRV(BR,BI,X,Y,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            ELSE IF (QI.EQ.4) THEN
               CALL C06FRU(BR,BI,X,Y,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            ELSE IF (QI.EQ.5) THEN
               CALL C06FRT(BR,BI,X,Y,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            ELSE IF (QI.EQ.6) THEN
               CALL C06FRS(BR,BI,X,Y,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            ELSE
               CALL C06FRR(BR,BI,X,Y,M*P,QI,R,TRIG((P-1)*QI*R),
     *                     TRIG(N+(P-1)*QI*R))
            END IF
         END IF
         INA = .NOT. INA
         P = P*QI
   20 CONTINUE
C
      FACTOR = 1.0D0/SQRT(DBLE(N))
      IF (INA) THEN
         DO 40 I = 0, M*N - 1
            X(I) = X(I)*FACTOR
            Y(I) = Y(I)*FACTOR
   40    CONTINUE
      ELSE
         DO 60 I = 0, M*N - 1
            X(I) = BR(I)*FACTOR
            Y(I) = BI(I)*FACTOR
   60    CONTINUE
      END IF
C
      RETURN
      END
