      SUBROUTINE F02AKY(N,LOW,IUPP,INTGER,HR,IHR,HI,IHI,VR,IVR,VI,IVI)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     COMPLEX ANALOGUE OF F01APF
C
C     .. Scalar Arguments ..
      INTEGER           IHI, IHR, IUPP, IVI, IVR, LOW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  HI(IHI,N), HR(IHR,N), VI(IVI,N), VR(IVR,N)
      INTEGER           INTGER(N)
C     .. Local Scalars ..
      INTEGER           I, I1, II, III, IUPP1, J, K, LOW1
C     .. Executable Statements ..
      DO 40 I = 1, N
         DO 20 J = 1, N
            VR(I,J) = 0.0D0
            VI(I,J) = 0.0D0
   20    CONTINUE
         VR(I,I) = 1.0D0
   40 CONTINUE
      IUPP1 = IUPP - 1
      LOW1 = LOW + 1
      DO 100 II = LOW1, IUPP1
         I = LOW1 + IUPP1 - II
         J = INTGER(I)
         I1 = I + 1
         DO 60 K = I1, IUPP
            III = I - 1
            VR(K,I) = HR(K,III)
            VI(K,I) = HI(K,III)
   60    CONTINUE
         IF (I.EQ.J) GO TO 100
         DO 80 K = I, IUPP
            VR(I,K) = VR(J,K)
            VI(I,K) = VI(J,K)
            VR(J,K) = 0.0D0
            VI(J,K) = 0.0D0
   80    CONTINUE
         VR(J,I) = 1.0D0
  100 CONTINUE
      RETURN
      END
