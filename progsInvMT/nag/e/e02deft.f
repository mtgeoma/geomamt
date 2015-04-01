      SUBROUTINE E02DEF(M,PX,PY,X,Y,LAMDA,MU,C,FF,WRK,IWRK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Derived from DASL routine B2VRE.
C
C     E02DEF. An algorithm for evaluating a bicubic polynomial
C     spline S(X,Y) from its B-spline representation at the M
C     points (X(I),Y(I)), I = 1, 2, ..., M.
C
C     Input Parameters:
C        M        The number of evaluation points.
C        PX       NXKNTS + 8,  where  NXKNTS  is number
C                    of interior  X-knots.
C        PY       NYKNTS + 8,  where  NYKNTS  is number
C                    of interior  Y-knots.
C        X        X-values.
C        Y        Y-values.
C        LAMDA    The X-knots.
C        MU       The Y-knots.
C        C        B-spline coefficients of  S.  First
C                    subscript relates to  X.
C
C     Output parameter:
C        FF       Values of spline.
C                 On exit, FF(I) contains the value of
C                 the spline evaluated at point (X(I),Y(I)),
C                 for I = 1,..,M.
C
C     Workspace parameters:
C        WRK      Real workspace of dimension at least (PY-4).
C        IWRK     Integer workspace of dimension at least (PY-4).
C
C     Failure indicator parameter:
C        IFAIL    Failure indicator:
C                 1 -  PX     .LT. 8,    or
C                      PY     .LT. 8,    or
C                      M      .LT. 1.
C                 2 -  E02DFW  failure for  X  or  Y.
C                 3 -  At least one point (X(K),Y(K)) lies
C                      outside the rectangle defined by
C                      LAMDA(4), LAMDA(PX-3), MU(4) and MU(PY-3).
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E02DEF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M, PX, PY
C     .. Array Arguments ..
      DOUBLE PRECISION  C((PX-4)*(PY-4)), FF(M), LAMDA(PX), MU(PY),
     *                  WRK(PY-4), X(M), Y(M)
      INTEGER           IWRK(PY-4)
C     .. Local Scalars ..
      DOUBLE PRECISION  XMAX, XMIN, YMAX, YMIN
      INTEGER           I, IERR, JERROR, K, NREC, NXKNTS, NYKNTS
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          E02DEZ, E02DFW
C     .. Executable Statements ..
      IERR = 0
      NREC = 0
      IF (M.LT.1 .OR. PX.LT.8 .OR. PY.LT.8) THEN
         IERR = 1
         NREC = 2
         WRITE (REC,FMT=99999) M, PX, PY
      ELSE
C        NXKNTS and NYKNTS are the number of interior knots in the
C        X and Y directions respectively.
         NXKNTS = PX - 8
         NYKNTS = PY - 8
C        Check X and Y knot sets.
         XMIN = LAMDA(4)
         XMAX = LAMDA(PX-3)
         YMIN = MU(4)
         YMAX = MU(PY-3)
         CALL E02DFW(PX,XMIN,XMAX,LAMDA,PX,JERROR)
         IF (JERROR.NE.0) THEN
            IERR = 2
            NREC = 1
            WRITE (REC,FMT=99998) 'LAMDA'
         ELSE
            CALL E02DFW(PY,YMIN,YMAX,MU,PY,JERROR)
            IF (JERROR.NE.0) THEN
               IERR = 2
               NREC = 1
               WRITE (REC,FMT=99998) 'MU'
            ELSE
C              Check that all points lie inside the spline domain.
               K = 0
               DO 20 I = M, 1, -1
                  IF (X(I).LT.LAMDA(4) .OR. X(I).GT.LAMDA(PX-3)
     *                .OR. Y(I).LT.MU(4) .OR. Y(I).GT.MU(PY-3)) THEN
                     K = I
                  END IF
   20          CONTINUE
               IF (K.GT.0) THEN
                  IERR = 3
                  NREC = 3
                  WRITE (REC,FMT=99997) K, X(K), Y(K)
               ELSE
C                 Evaluate S(X,Y) at all specified points.
                  CALL E02DEZ(M,NXKNTS,XMIN,XMAX,LAMDA,PX,NYKNTS,YMIN,
     *                        YMAX,MU,PY,C,PY-4,1,(PY-4)*(PX-4),X,Y,FF,
     *                        WRK,IWRK)
               END IF
            END IF
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, either M.lt.1, PX.lt.8 or PY.lt.8:',/4X,
     *       'M =',I16,', PX =',I16,', PY =',I16,'.')
99998 FORMAT (1X,'** On entry, the knots in ',A,' are not in non-decre',
     *       'asing order.')
99997 FORMAT (1X,'** On entry, point (X(K),Y(K)) lies outside the rect',
     *       'angle bounded by',/4X,'LAMDA(4), LAMDA(PX-3), MU(4), MU(',
     *       'PY-3):',/4X,'K =',I8,', X(K) = ',1P,D13.5,', Y(K) = ',
     *       D13.5,' .')
      END
