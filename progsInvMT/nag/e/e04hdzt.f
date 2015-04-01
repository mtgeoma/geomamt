      SUBROUTINE E04HDZ(N,LH,HESL,HESD,P,W,Y)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     FORMS THE PRODUCT Y = PT*H*P, WHERE P IS A VECTOR AND H IS AN
C     N*N SYMMETRIC MATRIX WITH ITS LOWER TRIANGLE STORED ROW BY ROW
C     IN HESL AND ITS DIAGONAL STORED IN HESD.
C
C     PHILIP E. GILL, ENID M. R. LONG, WALTER MURRAY AND SUSAN
C     PICKEN, D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND.
C     NOVEMBER 1976
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  Y
      INTEGER           LH, N
C     .. Array Arguments ..
      DOUBLE PRECISION  HESD(N), HESL(LH), P(N), W(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM, WI
      INTEGER           I, IM1, IP, IP1, IQ, K, L
C     .. Executable Statements ..
      L = 1
      DO 60 I = 1, N
         WI = HESD(I)*P(I)
         IF (I.EQ.1) GO TO 40
         IM1 = I - 1
         DO 20 K = 1, IM1
            WI = WI + HESL(L)*P(K)
            L = L + 1
   20    CONTINUE
   40    W(I) = WI
   60 CONTINUE
      L = 0
      SUM = 0.0D+0
      DO 120 I = 1, N
         L = L + I
         IP = I
         IQ = L
         WI = W(I)
         IF (I.EQ.N) GO TO 100
         IP1 = I + 1
         DO 80 K = IP1, N
            WI = WI + HESL(IQ)*P(K)
            IQ = IQ + IP
            IP = IP + 1
   80    CONTINUE
  100    SUM = SUM + WI*P(I)
  120 CONTINUE
      Y = SUM
      RETURN
C
C     END OF E04HDZ (PTHP)
C
      END
