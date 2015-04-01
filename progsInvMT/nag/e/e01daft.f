      SUBROUTINE E01DAF(MX,MY,X,Y,F,PX,PY,LAMDA,MU,C,WRK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     Derived from DASL routine B2IRE.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E01DAF')
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           IFAIL, MX, MY, PX, PY
C     .. Array Arguments ..
      DOUBLE PRECISION  C(MX*MY), F(MX*MY), LAMDA(MX+4), MU(MY+4),
     *                  WRK((MX+6)*(MY+6)), X(MX), Y(MY)
C     .. Local Scalars ..
      INTEGER           IDUM1, IDUM2, IE, IER2, IERR, IWGHT, IWRK,
     *                  IXKNOT, IXROW, IXU, IYKNOT, IYROW, IYU, JERROR,
     *                  LWRK, MDIST, NE, NREC, NXKNTS, NXU, NYKNTS, NYU
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          E01DAX, E01DAY, E01DAZ
C     .. Executable Statements ..
      IERR = 0
      NREC = 0
      IF (MX.LT.4) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99999) MX
      ELSE IF (MY.LT.4) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99998) MY
      ELSE
C        Numbers of interior X- and Y-knots.
         NXKNTS = MX - 4
         NYKNTS = MY - 4
C        Number of elements within intermediate coefficient matrix E.
         NE = MX*MY
C        Numbers of elements within bands of band upper triangular
C        matrices XUFCTR and YUFCTR.
         NXU = 2*(MX+NXKNTS+1)
         NYU = 2*(MY+NYKNTS+1)
C        Subdivide workspace.
         IE = 1
         IXU = IE + NE
         IYU = IXU + NXU
         IXKNOT = IYU + NYU
         IYKNOT = IXKNOT + 8
         IXROW = IYKNOT + 8
         IYROW = IXROW + MX
         IWRK = IYROW + MY
C        Check endpoints and X- and Y-values.
         MDIST = 0
         IWGHT = 1
         WRK(IWGHT) = ONE
         CALL E01DAX(X(1),X(MX),MX,X,WRK(IWGHT),0,1,IDUM1,IDUM2,MDIST,
     *               JERROR)
         IF (JERROR.EQ.0 .AND. MDIST.EQ.MX) THEN
            CALL E01DAX(Y(1),Y(MY),MY,Y,WRK(IWGHT),0,1,IDUM1,IDUM2,
     *                  MDIST,JERROR)
            IF (JERROR.EQ.0 .AND. MDIST.EQ.MY) THEN
C              Place X- and Y-knotlines.
               CALL E01DAY(4,NXKNTS,MX,X,WRK(IWGHT),0,1,1,MX,LAMDA(5),
     *                     MX)
               CALL E01DAY(4,NYKNTS,MY,Y,WRK(IWGHT),0,1,1,MY,MU(5),MY)
C              Construct spline interpolant.
               LWRK = (MX+6)*(MY+6)
               CALL E01DAZ(4,NXKNTS,X(1),X(MX),4,NYKNTS,Y(1),Y(MY),MX,X,
     *                     MY,Y,F,MY,1,MY*MX,LAMDA(5),MX,MU(5),MY,C,MY,
     *                     1,MY*MX,WRK(IXU),NXU,WRK(IYU),NYU,WRK(IE),1,
     *                     MX,NE,WRK(IXKNOT),8,WRK(IYKNOT),8,WRK(IXROW),
     *                     MX,WRK(IYROW),MY,WRK(IWRK),LWRK-IWRK+1,IER2)
C              Indicate whether system is numerically singular.
               IF (IER2.NE.0) THEN
                  IERR = 3
                  NREC = 2
                  WRITE (REC,FMT=99996)
               END IF
            ELSE
               IERR = 2
               NREC = 1
               WRITE (REC,FMT=99997)
            END IF
         ELSE
            IERR = 2
            NREC = 1
            WRITE (REC,FMT=99997)
         END IF
         IF (IERR.EQ.0) THEN
C           Insert exterior knots into LAMDA and MU.
            LAMDA(1) = X(1)
            LAMDA(2) = X(1)
            LAMDA(3) = X(1)
            LAMDA(4) = X(1)
            LAMDA(MX+1) = X(MX)
            LAMDA(MX+2) = X(MX)
            LAMDA(MX+3) = X(MX)
            LAMDA(MX+4) = X(MX)
            MU(1) = Y(1)
            MU(2) = Y(1)
            MU(3) = Y(1)
            MU(4) = Y(1)
            MU(MY+1) = Y(MY)
            MU(MY+2) = Y(MY)
            MU(MY+3) = Y(MY)
            MU(MY+4) = Y(MY)
C           Return the number of knots in X and Y directions, PX and PY.
            PX = MX + 4
            PY = MY + 4
         END IF
      END IF
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (1X,'** On entry, MX.lt.4: MX =',I16)
99998 FORMAT (1X,'** On entry, MY.lt.4: MY =',I16)
99997 FORMAT (1X,'** On entry, the X or the Y mesh points are not in s',
     *       'trictly ascending order.')
99996 FORMAT (1X,'** An intermediate set of linear equations is singul',
     *       'ar',/4X,'- the data is too ill-conditioned to compute B-',
     *       'spline coefficients.')
      END
