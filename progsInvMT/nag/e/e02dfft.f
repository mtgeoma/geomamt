      SUBROUTINE E02DFF(MX,MY,PX,PY,X,Y,LAMDA,MU,C,FF,WRK,LWRK,IWRK,
     *                  LIWRK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Derived from DASL routine B2VRE.
C
C     E02DFF.  An algorithm for evaluating a bivariate
C     polynomial spline  S(X,Y)  from its B-spline
C     representation at all vertices of a rectangular mesh.
C
C     The routine is designed to be fast, but this is at the
C     expense of some workspace.  If speed is less important
C     than space,  E02DFF  may be called several times with
C     one meshline of values, for example, at each call.
C
C     Input Parameters:
C        PX       NXKNTS + 8,  where  NXKNTS  is number
C                    of interior  X-knots.
C        LAMDA    The X-knots.
C        PY       NYKNTS + 8,  where  NYKNTS  is number
C                    of interior  Y-knots.
C        MU       The Y-knots.
C        C        B-spline coefficients of  S.  First
C                    subscript relates to  X.
C        MX       The number of  X-meshlines.
C        X        X-meshline values.
C        MY       The number of  Y-meshlines.
C        Y        The Y-meshline values.
C
C     Output parameters:
C        FF       Values of spline.
C                 On exit, FF((I-1)*MY+J) contains the value of
C                 the spline evaluated at point (X(I),Y(J)),
C                 for I = 1,..,MX, J = 1,..,MY.
C
C     Workspace (and associated dimension) parameters:
C        WRK      Real workspace.
C        LWRK     Dimension of  WRK.  .ge. NWRK =
C                    min(NWRK1, NWRK2),  where
C                    NWRK1 = MX*4 + PX,
C                    NWRK2 = MY*4 + PY.
C        IWRK     Integer workspace.
C        LIWRK    Dimension of  IWRK.  .ge. NIWRK =
C                    MX + PX - 4  (if NWRK1 .le. NWRK2),
C                    MY + PY - 4  (otherwise).
C
C     Failure indicator parameter:
C        IFAIL    Failure indicator:
C                 1 -  PX     .LT. 8,    or
C                      PY     .LT. 8,    or
C                      MX     .LT. 0,    or
C                      MY     .LT. 0,    or
C                      LWRK   .LT. NWRK  or
C                      LIWRK  .LT. NIWRK
C                 2 -  E02DFW  failure for  X  or  Y.
C                 3 -  LAMDA(4) .le. X(1) .lt. X(2) .lt. ...
C                               .lt. X(MX) .le. LAMDA(PX-3)    or
C                      MU(4) .le. Y(1) .lt. Y(2) .lt. ...
C                            .lt. Y(MY) .le. MU(PY-3)
C                      violated.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E02DFF')
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LIWRK, LWRK, MX, MY, PX, PY
C     .. Array Arguments ..
      DOUBLE PRECISION  C((PX-4)*(PY-4)), FF(MX*MY), LAMDA(PX), MU(PY),
     *                  WRK(LWRK), X(MX), Y(MY)
      INTEGER           IWRK(LIWRK)
C     .. Local Scalars ..
      DOUBLE PRECISION  XMAX, XMIN, YMAX, YMIN
      INTEGER           IDUMMY, IERR, JERROR, NIWRK, NREC, NWRK, NWRK1,
     *                  NWRK2, NXDIST, NXKNTS, NYDIST, NYKNTS
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          E01DAX, E02DFW, E02DFZ
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
      IERR = 0
      NREC = 0
      IF (PX.LT.8 .OR. PY.LT.8 .OR. MX.LT.1 .OR. MY.LT.1) THEN
         IERR = 1
         NREC = 2
         WRITE (REC,FMT=99999) PX, PY, MX, MY
      ELSE
C        NXKNTS and NYKNTS are the number of interior knots in the
C        X and Y directions respectively.
         NXKNTS = PX - 8
         NYKNTS = PY - 8
C        Real and integer workspace requirements.
         NWRK1 = MX*4 + PX
         NWRK2 = MY*4 + PY
         NWRK = MIN(NWRK1,NWRK2)
         IF (NWRK1.LE.NWRK2) THEN
            NIWRK = MX + PX - 4
         ELSE
            NIWRK = MY + PY - 4
         END IF
C        Basic checks on workspace requirements.
         IF (LWRK.LT.NWRK) THEN
            IERR = 2
            NREC = 1
            WRITE (REC,FMT=99998) NWRK, LWRK
         ELSE IF (LIWRK.LT.NIWRK) THEN
            IERR = 2
            NREC = 1
            WRITE (REC,FMT=99997) NIWRK, LIWRK
         ELSE
C           Check X and Y knot sets.
            XMIN = LAMDA(4)
            XMAX = LAMDA(PX-3)
            YMIN = MU(4)
            YMAX = MU(PY-3)
            CALL E02DFW(PX,XMIN,XMAX,LAMDA,PX,JERROR)
            IF (JERROR.NE.0) THEN
               IERR = 3
               NREC = 1
               WRITE (REC,FMT=99996) 'LAMDA'
            ELSE
               CALL E02DFW(PY,YMIN,YMAX,MU,PY,JERROR)
               IF (JERROR.NE.0) THEN
                  IERR = 3
                  NREC = 1
                  WRITE (REC,FMT=99996) 'MU'
               ELSE
C                 Check that the  X-values  are strictly increasing ...
                  WRK(1) = ONE
                  NXDIST = 0
                  CALL E01DAX(XMIN,XMAX,MX,X,WRK(1),0,1,IDUMMY,IDUMMY,
     *                        NXDIST,JERROR)
                  IF (NXDIST.NE.MX .OR. JERROR.NE.0) THEN
                     IERR = 4
                     NREC = 2
                     WRITE (REC,FMT=99995)
                  ELSE
C                    ... and then a similar check for Y.
                     NYDIST = 0
                     CALL E01DAX(YMIN,YMAX,MY,Y,WRK(1),0,1,IDUMMY,
     *                           IDUMMY,NYDIST,JERROR)
                     IF (NYDIST.NE.MY .OR. JERROR.NE.0) THEN
                        IERR = 4
                        NREC = 2
                        WRITE (REC,FMT=99994)
                     ELSE
C                       Evaluate  S(X, Y)  at all specified points.
                        IF (NWRK1.LE.NWRK2) THEN
                           CALL E02DFZ(4,NYKNTS,YMIN,YMAX,MU,PY,4,
     *                                 NXKNTS,XMIN,XMAX,LAMDA,PX,C,1,
     *                                 PY-4,(PX-4)*(PY-4),MY,Y,MX,X,FF,
     *                                 1,MY,MY*MX,WRK,LWRK,IWRK,LIWRK)
                        ELSE
                           CALL E02DFZ(4,NXKNTS,XMIN,XMAX,LAMDA,PX,4,
     *                                 NYKNTS,YMIN,YMAX,MU,PY,C,PY-4,1,
     *                                 (PY-4)*(PX-4),MX,X,MY,Y,FF,MY,1,
     *                                 MX*MY,WRK,LWRK,IWRK,LIWRK)
                        END IF
                     END IF
                  END IF
               END IF
            END IF
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, either PX.lt.8, PY.lt.8, MX.lt.1 or MY.',
     *       'lt.1:',/4X,'PX =',I13,', PY =',I13,', MX =',I13,', MY =',
     *       I13,'.')
99998 FORMAT (1X,'** On entry, LWRK.lt.',I8,': LWRK =',I16)
99997 FORMAT (1X,'** On entry, LIWRK.lt.',I8,': LIWRK =',I16)
99996 FORMAT (1X,'** On entry, the knots in ',A,' are not in non-decre',
     *       'asing order.')
99995 FORMAT (1X,'** On entry, the restriction  LAMDA(4).le.X(1).lt. .',
     *       '.. .lt.X(MX).le.LAMDA(PX-3)',/4X,'is violated.')
99994 FORMAT (1X,'** On entry, the restriction  MU(4).le.Y(1).lt. ... ',
     *       '.lt.Y(MY).le.MU(PY-3)',/4X,'is violated.')
      END
