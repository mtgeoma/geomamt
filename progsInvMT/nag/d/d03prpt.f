      SUBROUTINE D03PRP(XP,UP,IPTS,X,M,U,NPTS,NPDE,IFAIL,C,IXFIX,NFIX)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C***********************************************************************
C     Sprint cubic spline interpolation routine, driver for piece-poly
C
C     N.B. This routine assumes that the array ifix contains the indices
C        of the fixed points of the existing mesh
C     I.E.  X(IXFIX(I))  is a fixed point at which there may be a PDE
C        discontinuity.
C
C      INTDRV routine from SPRINT
C***********************************************************************
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IPTS, M, NFIX, NPDE, NPTS
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NPTS,3), U(NPDE,NPTS), UP(NPDE,IPTS), X(NPTS),
     *                  XP(IPTS)
      INTEGER           IXFIX(*)
C     .. Scalars in Common ..
      INTEGER           IDEV, ITRACE
C     .. Local Scalars ..
      INTEGER           I, IT, J, JJ, JK, K, KK, NKPTS, NXFXP1
C     .. External Subroutines ..
      EXTERNAL          D03PRT, D03PZW
C     .. Common blocks ..
      COMMON            /AD02NM/ITRACE, IDEV
C     .. Save statement ..
      SAVE              /AD02NM/
C     .. Executable Statements ..
      IF (NFIX.EQ.0) THEN
C         No material interfaces present -  use ordinary cubic spline.
         CALL D03PRT(XP,UP,IPTS,X,U,NPTS,NPDE,IFAIL,C)
      ELSE
C         Cubic spline routine must be called in piecewise fashion .
         NXFXP1 = NFIX + 1
         JJ = 1
         DO 60 I = 1, NXFXP1
            JK = 0
            IF (I.EQ.1) THEN
               K = 1
            ELSE
               K = IXFIX(I-1) + 0.1D0
            END IF
C            X(K) is the starting point for the next element
            IF (I.EQ.NXFXP1) THEN
               NKPTS = NPTS - K + 1
               KK = NPTS
            ELSE
               NKPTS = IXFIX(I) - K + 1.1D0
               KK = IXFIX(I) + 0.1D0
            END IF
   20       J = JJ + JK
            IF (J.GT.IPTS) THEN
               GO TO 40
            ELSE
               IF (XP(J).LE.X(KK)) THEN
                  JK = JK + 1
C               Counts the no. of new points in this element
                  GO TO 20
               END IF
            END IF
   40       J = J - 1
C           Can interpolate to find new solution values from
C           XP(JJ)th mesh point to XP(J)th mesh point.
C       **************  IF (ITRACE.GE.1)WRITE(IDEV,90)JJ,J,JK,K,NKPTS
C    90          FORMAT(' XP(',I3,') TO XP(',I3,' USING',I4,'POINTS'/
C       ****1********   ' INTERP FROM X(',I3,') USING ',I4,' POINTS')
            IF (JK.GT.0 .AND. NKPTS.GE.4) THEN
               CALL D03PRT(XP(JJ),UP(1,JJ),JK,X(K),U(1,K),NKPTS,NPDE,
     *                     IFAIL,C)
            ELSE IF (JK.GT.0) THEN
               IT = 1
               CALL D03PZW(XP(JJ),UP(1,JJ),JK,X(K),M,U(1,K),NKPTS,NPDE,
     *                     IT,IFAIL)
            END IF
            JJ = J
   60    CONTINUE
      END IF
      RETURN
      END
