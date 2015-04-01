      SUBROUTINE F02BDF(AR,IAR,AI,IAI,N,RLB,RUB,M,MM,RR,RI,VR,IVR,VI,
     *                  IVI,INTGER,C,BR,IBR,BI,IBI,U,V,LFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     W.PHILLIPS. OXFORD UNIVERSITY COMPUTING SERVICE. 1 JUN 1977
C
C     SELECTED EIGENVALUES AND EIGENVECTORS OF A COMPLEX
C     UNSYMMETRIC MATRIX
C
C     THIS ROUTINE REPLACES F02ALF.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02BDF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RLB, RUB
      INTEGER           IAI, IAR, IBI, IBR, IVI, IVR, LFAIL, M, MM, N
C     .. Array Arguments ..
      DOUBLE PRECISION  AI(IAI,N), AR(IAR,N), BI(IBI,N), BR(IBR,N),
     *                  RI(N), RR(N), U(N), V(N), VI(IVI,M), VR(IVR,M)
      INTEGER           INTGER(N)
      LOGICAL           C(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC, CT, DT, RMAX, RRI, RRR, SUM, TERM, XXXX,
     *                  ZERO
      INTEGER           I, ISAVE, J
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  A02ABF, X02AJF
      INTEGER           P01ABF
      EXTERNAL          A02ABF, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01AMF, F01ANF, F02ANF, F02BLF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Data statements ..
      DATA              ZERO/0.0D0/
C     .. Executable Statements ..
      ISAVE = LFAIL
      LFAIL = 1
      CALL F01AMF(N,1,N,AR,IAR,AI,IAI,INTGER)
      ACC = X02AJF()
      DO 40 I = 1, N
         DO 20 J = 1, N
            BR(I,J) = AR(I,J)
            BI(I,J) = AI(I,J)
   20    CONTINUE
   40 CONTINUE
      CALL F02ANF(N,ACC,BR,IBR,BI,IBI,RR,RI,LFAIL)
      IF (LFAIL.EQ.0) GO TO 60
      LFAIL = P01ABF(ISAVE,LFAIL,SRNAME,0,P01REC)
      RETURN
   60 MM = 0
      DO 80 I = 1, N
         C(I) = .FALSE.
   80 CONTINUE
      DO 100 I = 1, N
         IF (A02ABF(RR(I),RI(I)).LT.RLB) GO TO 100
         IF (A02ABF(RR(I),RI(I)).GT.RUB) GO TO 100
         MM = MM + 1
         C(I) = .TRUE.
  100 CONTINUE
      IF (MM.NE.0) GO TO 120
      LFAIL = P01ABF(ISAVE,2,SRNAME,0,P01REC)
      RETURN
  120 IF (MM.LE.M) GO TO 140
      LFAIL = P01ABF(ISAVE,3,SRNAME,0,P01REC)
      RETURN
  140 CALL F02BLF(N,MM,AR,IAR,AI,IAI,RI,C,RR,VR,IVR,VI,IVI,BR,IBR,BI,
     *            IBI,U,V,LFAIL)
      CALL F01ANF(1,N,MM,AR,IAR,AI,IAI,INTGER,VR,IVR,VI,IVI,N)
      J = 0
      DO 160 I = 1, N
         IF ( .NOT. C(I)) GO TO 160
         J = J + 1
         IF (J.EQ.I) GO TO 160
         RRR = RR(J)
         RRI = RI(J)
         RR(J) = RR(I)
         RI(J) = RI(I)
         RR(I) = RRR
         RI(I) = RRI
  160 CONTINUE
      DO 260 I = 1, MM
         SUM = ZERO
         RMAX = ZERO
         DO 180 J = 1, N
            IF (ABS(VR(J,I)).GT.RMAX) RMAX = ABS(VR(J,I))
            IF (ABS(VI(J,I)).GT.RMAX) RMAX = ABS(VI(J,I))
  180    CONTINUE
         DO 200 J = 1, N
            VR(J,I) = VR(J,I)/RMAX
            VI(J,I) = VI(J,I)/RMAX
  200    CONTINUE
         RMAX = ZERO
         DO 220 J = 1, N
            TERM = VR(J,I)**2 + VI(J,I)**2
            SUM = SUM + TERM
            IF (TERM.LE.RMAX) GO TO 220
            RMAX = TERM
            CT = VR(J,I)
            DT = -VI(J,I)
  220    CONTINUE
         IF (SUM.EQ.ZERO) GO TO 260
         SUM = SQRT(SUM*RMAX)
         DO 240 J = 1, N
            TERM = VR(J,I)
            VR(J,I) = (VR(J,I)*CT-VI(J,I)*DT)/SUM
            VI(J,I) = (DT*TERM+CT*VI(J,I))/SUM
  240    CONTINUE
  260 CONTINUE
      RETURN
      END
