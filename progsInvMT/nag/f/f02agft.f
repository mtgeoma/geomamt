      SUBROUTINE F02AGF(A,IA,N,RR,RI,VR,IVR,VI,IVI,INTGER,IFAIL)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     MARK 14A REVISED. IER-685 (DEC 1989).
C
C     EIGENVALUES AND EIGENVECTORS OF REAL UNSYMMETRIC MATRIX
C     1ST AUGUST 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02AGF')
C     .. Scalar Arguments ..
      INTEGER           IA, IFAIL, IVI, IVR, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), RI(N), RR(N), VI(IVI,N), VR(IVR,N)
      INTEGER           INTGER(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  C, D, MACHEP, MAX, SUM, TERM
      INTEGER           I, IB, ISAVE, J, K, L
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF, X02BHF
      EXTERNAL          X02AJF, P01ABF, X02BHF
C     .. External Subroutines ..
      EXTERNAL          F01AKF, F01APF, F01ATF, F01AUF, F02AQF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL = 1
      MACHEP = X02AJF()
      IB = X02BHF()
      CALL F01ATF(N,IB,A,IA,K,L,RR)
      CALL F01AKF(N,K,L,A,IA,INTGER)
      CALL F01APF(N,K,L,INTGER,A,IA,VR,IVR)
      CALL F01AUF(N,K,L,N,RR,VR,IVR)
      CALL F02AQF(N,1,N,MACHEP,A,IA,VR,IVR,RR,RI,INTGER,IFAIL)
      IF (IFAIL.EQ.0) GO TO 20
      IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
   20 DO 140 I = 1, N
         IF (RI(I).EQ.0.0D0) GO TO 60
         IF (RI(I).GT.0.0D0) GO TO 100
         DO 40 J = 1, N
            VR(J,I) = VR(J,I-1)
            VI(J,I) = -VI(J,I-1)
   40    CONTINUE
         GO TO 140
   60    DO 80 J = 1, N
            VI(J,I) = 0.0D0
   80    CONTINUE
         GO TO 140
  100    DO 120 J = 1, N
            VI(J,I) = VR(J,I+1)
  120    CONTINUE
  140 CONTINUE
      DO 280 I = 1, N
         SUM = 0.0D0
         MAX = 0.0D0
         DO 180 J = 1, N
            IF (ABS(VR(J,I)).LE.MAX) GO TO 160
            MAX = ABS(VR(J,I))
  160       IF (ABS(VI(J,I)).LE.MAX) GO TO 180
            MAX = ABS(VI(J,I))
  180    CONTINUE
         DO 200 J = 1, N
            VR(J,I) = VR(J,I)/MAX
            VI(J,I) = VI(J,I)/MAX
  200    CONTINUE
         MAX = 0.0D0
         DO 240 J = 1, N
            TERM = VR(J,I)**2 + VI(J,I)**2
            SUM = SUM + TERM
            IF (TERM.LE.MAX) GO TO 220
            MAX = TERM
            C = VR(J,I)
            D = -VI(J,I)
  220       CONTINUE
  240    CONTINUE
         SUM = SUM*(C**2+D**2)
         SUM = SQRT(SUM)
         DO 260 J = 1, N
            TERM = VR(J,I)
            VR(J,I) = (VR(J,I)*C-VI(J,I)*D)/SUM
            VI(J,I) = (D*TERM+C*VI(J,I))/SUM
  260    CONTINUE
  280 CONTINUE
      RETURN
      END
