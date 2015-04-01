      SUBROUTINE D03PZW(XP,UP,IPTS,X,M,U,NPTS,NPDE,ITYPE,IFAIL1)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ---------------------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER           IFAIL1, IPTS, ITYPE, M, NPDE, NPTS
C     .. Array Arguments ..
      DOUBLE PRECISION  U(NPDE,*), UP(NPDE,IPTS,ITYPE), X(NPTS),
     *                  XP(IPTS)
C     .. Scalars in Common ..
      INTEGER           NIJ
C     .. Local Scalars ..
      DOUBLE PRECISION  DL, DPL, DPR, DR, H, HH, HSUM, PL, PR, UL, UM,
     *                  UR, UROUND, XHAT, XL, XM, XR, XSUM
      INTEGER           I, J, JI, NC
      LOGICAL           SINGLR
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. Intrinsic Functions ..
      INTRINSIC         LOG
C     .. Common blocks ..
      COMMON            /AD03PZ/NIJ
C     .. Executable Statements ..
C
      UROUND = X02AJF()
      NC = M
      IF (X(1).LT.UROUND .AND. NC.GT.0) THEN
         SINGLR = .TRUE.
      ELSE
         SINGLR = .FALSE.
      END IF
C
C     ... A better check would to test when the problem is run ...
C
      DO 100 JI = 1, IPTS
         IFAIL1 = 0
C
C        ... Search the grid points which bound XP(JI) ...
C
         IF (XP(JI).LT.X(1)-UROUND) GO TO 40
         DO 20 I = 2, NPTS
            IF (XP(JI).LE.X(I)) GO TO 60
   20    CONTINUE
         IF (XP(JI).GT.X(NPTS)+UROUND) GO TO 40
         I = NPTS
         GO TO 60
C
C        ... XP(JI) is outside the mesh ...
C
   40    CONTINUE
         IFAIL1 = 3
         NIJ = JI
         GO TO 120
   60    CONTINUE
C
C        ... The point XP(JI) is bounded by the i-1th ...
C        ... and the ith grid pts ...
C
         XL = X(I-1)
         XR = X(I)
         H = XR - XL
         XHAT = XP(JI)
         XSUM = XR + XL
C
C        ... Non-polar case ...
C
         IF (NC.EQ.0) THEN
            PR = (XHAT-XL)/H
            DPR = 1.D0/H
            DPL = -DPR
         ELSE IF (SINGLR) THEN
            PR = (XHAT**2-XL**2)/(XSUM*H)
            DPR = 2*XHAT/(XSUM*H)
            DPL = -DPR
         ELSE IF (NC.EQ.1) THEN
C
C           ... Cylindrical polar co-ordinates ...
C
            DPR = 1.D0/LOG(XR/XL)
            PR = LOG(XHAT/XL)*DPR
            DPR = DPR/XHAT
            DPL = -DPR
         ELSE IF (NC.EQ.2) THEN
            PR = XR*(1.D0-XL/XHAT)/H
            DPR = XL*XR/(XHAT**2*H)
            DPL = DPR
         END IF
         PL = 1.0D0 - PR
         DO 80 J = 1, NPDE
C
C           ... Apply linear interpolation ...
C
            UP(J,JI,1) = U(J,I)*PR + U(J,I-1)*PL
            IF (ITYPE.GE.2) THEN
C
C              ... Compute ist space deriv at XP(JI) ...
C
               IF (XP(JI).LT.XR .OR. I.EQ.NPTS) THEN
                  UP(J,JI,2) = U(J,I)*DPR + U(J,I-1)*DPL
               ELSE
C
C                 ... Interpolation at an interior mesh point
C
                  UL = U(J,I-1)
                  UM = U(J,I)
                  UR = U(J,I+1)
                  XM = XR
                  XR = X(I+1)
                  HH = XR - XM
                  HSUM = HH + H
C
C                 ... Linear interpolation ...
C
                  IF (NC.EQ.0) THEN
                     UP(J,JI,2) = (HH*(UM-UL)/H+H*(UR-UM)/HH)/HSUM
C
C                    ... Quadratic interp for singular problem ...
C
                  ELSE IF (SINGLR) THEN
                     UP(J,JI,2) = 2.D0*XM*(HH*(UM-UL)/(H*(XL+XM))
     *                            +H*(UR-UM)/(HH*(XR+XM)))/HSUM
C
C                    ... Use log form of basis functions ...
C
                  ELSE IF (NC.EQ.1) THEN
                     DL = LOG(XM/XL)
                     DR = LOG(XR/XM)
                     UP(J,JI,2) = (HH*(UM-UL)/DL+H*(UR-UM)/DR)/(XM*HSUM)
                  ELSE IF (NC.EQ.2) THEN
                     UP(J,JI,2) = (HH*XL*XM*(UM-UL)/H+H*XR*XM*(UR-UM)
     *                            /HH)/(XM**2*HSUM)
                  END IF
               END IF
            END IF
   80    CONTINUE
  100 CONTINUE
  120 CONTINUE
C
      RETURN
      END
