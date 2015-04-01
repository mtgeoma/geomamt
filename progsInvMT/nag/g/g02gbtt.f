      SUBROUTINE G02GBT(MEAN,N,M,X,LDX,ISX,IP,Q,LDQ,SVD,IRANK,WWT,H,WK)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Scalar Arguments ..
      INTEGER           IP, IRANK, LDQ, LDX, M, N
      LOGICAL           SVD
      CHARACTER         MEAN
C     .. Array Arguments ..
      DOUBLE PRECISION  H(N), Q(LDQ,*), WK(2*IP), WWT(N), X(LDX,M)
      INTEGER           ISX(M)
C     .. Local Scalars ..
      INTEGER           I, IM, J, K
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DTRSV
C     .. Executable Statements ..
      IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') THEN
         IM = 1
      ELSE
         IM = 0
      END IF
      IF (SVD) THEN
         WK(1) = 1.0D0
         DO 40 I = 1, N
            K = IM
            DO 20 J = 1, M
               IF (ISX(J).GT.0) THEN
                  K = K + 1
                  WK(K) = X(I,J)
               END IF
   20       CONTINUE
            CALL DGEMV('N',IRANK,IP,1.0D0,Q,LDQ,WK,1,0.0D0,WK(IP+1),1)
            H(I) = DDOT(IRANK,WK(IP+1),1,WK(IP+1),1)*WWT(I)*WWT(I)
   40    CONTINUE
      ELSE
         DO 80 I = 1, N
            WK(1) = 1.0D0
            K = IM
            DO 60 J = 1, M
               IF (ISX(J).GT.0) THEN
                  K = K + 1
                  WK(K) = X(I,J)
               END IF
   60       CONTINUE
            CALL DTRSV('U','T','N',IP,Q,LDQ,WK,1)
            H(I) = DDOT(IP,WK,1,WK,1)*WWT(I)*WWT(I)
   80    CONTINUE
      END IF
      RETURN
      END
