      SUBROUTINE F02AKF(AR,IAR,AI,IAI,N,WR,WI,VR,IVR,VI,IVI,INTGER,
     *                  IFAIL)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     MARK 14A REVISED. IER-686 (DEC 1989).
C
C     EIGENVALUES AND EIGENVECTORS OF A COMPLEX MATRIX
C     1ST AUGUST 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02AKF')
C     .. Scalar Arguments ..
      INTEGER           IAI, IAR, IFAIL, IVI, IVR, N
C     .. Array Arguments ..
      DOUBLE PRECISION  AI(IAI,N), AR(IAR,N), VI(IVI,N), VR(IVR,N),
     *                  WI(N), WR(N)
      INTEGER           INTGER(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACHEPS, C, D, MAX, SUM, TERM
      INTEGER           I, IB, IERR, J, K, L
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF, X02BHF
      EXTERNAL          X02AJF, P01ABF, X02BHF
C     .. External Subroutines ..
      EXTERNAL          F01AMF, F01AVF, F01AWF, F02AKY, F02AKZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      IB = X02BHF()
      CALL F01AVF(N,IB,AR,IAR,AI,IAI,K,L,WR)
      CALL F01AMF(N,K,L,AR,IAR,AI,IAI,INTGER)
      ACHEPS = X02AJF()
      CALL F02AKY(N,K,L,INTGER,AR,IAR,AI,IAI,VR,IVR,VI,IVI)
      CALL F01AWF(N,K,L,N,WR,VR,IVR,VI,IVI)
      CALL F02AKZ(N,1,N,ACHEPS,AR,IAR,AI,IAI,WR,WI,VR,IVR,VI,IVI,IERR)
      IF (IERR.EQ.0) GO TO 20
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
   20 DO 160 I = 1, N
         SUM = 0.0D0
         MAX = 0.0D0
         DO 60 J = 1, N
            IF (ABS(VR(J,I)).LE.MAX) GO TO 40
            MAX = ABS(VR(J,I))
   40       IF (ABS(VI(J,I)).LE.MAX) GO TO 60
            MAX = ABS(VI(J,I))
   60    CONTINUE
         DO 80 J = 1, N
            VR(J,I) = VR(J,I)/MAX
            VI(J,I) = VI(J,I)/MAX
   80    CONTINUE
         MAX = 0.0D0
         DO 120 J = 1, N
            TERM = VR(J,I)**2 + VI(J,I)**2
            SUM = SUM + TERM
            IF (TERM.LE.MAX) GO TO 100
            MAX = TERM
            C = VR(J,I)
            D = -VI(J,I)
  100       CONTINUE
  120    CONTINUE
         SUM = SUM*(C**2+D**2)
         SUM = SQRT(SUM)
         DO 140 J = 1, N
            TERM = VR(J,I)
            VR(J,I) = (VR(J,I)*C-VI(J,I)*D)/SUM
            VI(J,I) = (D*TERM+C*VI(J,I))/SUM
  140    CONTINUE
  160 CONTINUE
      IFAIL = 0
      RETURN
      END
