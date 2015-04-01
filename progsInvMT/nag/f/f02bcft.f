      SUBROUTINE F02BCF(A,IA,N,RLB,RUB,M,MM,RR,RI,VR,IVR,VI,IVI,INTGER,
     *                  ICNT,C,B,IB,U,V,LFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     W.PHILLIPS. OXFORD UNIVERSITY COMPUTING SERVICE. 1 JUN 1977
C
C     SELECTED EIGENVALUES AND EIGENVECTORS  OF A REAL UNSYMMETRIC
C     MATRIX
C
C     THIS ROUTINE REPLACES F02AHF.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02BCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RLB, RUB
      INTEGER           IA, IB, IVI, IVR, LFAIL, M, MM, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), B(IB,N), RI(N), RR(N), U(N), V(N),
     *                  VI(IVI,M), VR(IVR,M)
      INTEGER           ICNT(N), INTGER(N)
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
      EXTERNAL          F01AKF, F01ALF, F02APF, F02BKF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Data statements ..
      DATA              ZERO/0.0D0/
C     .. Executable Statements ..
      ISAVE = LFAIL
      LFAIL = 1
      CALL F01AKF(N,1,N,A,IA,INTGER)
      ACC = X02AJF()
      DO 40 I = 1, N
         DO 20 J = 1, N
            B(I,J) = A(I,J)
   20    CONTINUE
   40 CONTINUE
      CALL F02APF(N,ACC,B,IB,RR,RI,ICNT,LFAIL)
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
  140 CALL F02BKF(N,MM,A,IA,RI,C,RR,VR,IVR,B,IB,U,V,LFAIL)
      CALL F01ALF(1,N,MM,A,IA,INTGER,VR,IVR,N)
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
      DO 280 I = 1, MM
         IF (RI(I).EQ.ZERO) GO TO 200
         IF (RI(I).GT.ZERO) GO TO 240
         DO 180 J = 1, N
            VR(J,I) = VR(J,I-1)
            VI(J,I) = -VI(J,I-1)
  180    CONTINUE
         GO TO 280
  200    DO 220 J = 1, N
            VI(J,I) = ZERO
  220    CONTINUE
         GO TO 280
  240    DO 260 J = 1, N
            VI(J,I) = VR(J,I+1)
  260    CONTINUE
  280 CONTINUE
      DO 380 I = 1, MM
         SUM = ZERO
         RMAX = ZERO
         DO 300 J = 1, N
            IF (ABS(VR(J,I)).GT.RMAX) RMAX = ABS(VR(J,I))
            IF (ABS(VI(J,I)).GT.RMAX) RMAX = ABS(VI(J,I))
  300    CONTINUE
         DO 320 J = 1, N
            VR(J,I) = VR(J,I)/RMAX
            VI(J,I) = VI(J,I)/RMAX
  320    CONTINUE
         RMAX = ZERO
         DO 340 J = 1, N
            TERM = VR(J,I)**2 + VI(J,I)**2
            SUM = SUM + TERM
            IF (TERM.LE.RMAX) GO TO 340
            RMAX = TERM
            CT = VR(J,I)
            DT = -VI(J,I)
  340    CONTINUE
         IF (SUM.EQ.ZERO) GO TO 380
         SUM = SQRT(SUM*RMAX)
         DO 360 J = 1, N
            TERM = VR(J,I)
            VR(J,I) = (VR(J,I)*CT-VI(J,I)*DT)/SUM
            VI(J,I) = (DT*TERM+CT*VI(J,I))/SUM
  360    CONTINUE
  380 CONTINUE
      RETURN
      END
