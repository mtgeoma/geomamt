      SUBROUTINE D03EEF(XMIN,XMAX,YMIN,YMAX,PDEF,BNDY,NGX,NGY,LDA,A,RHS,
     *                  SCHEME,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D03EEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XMAX, XMIN, YMAX, YMIN
      INTEGER           IFAIL, LDA, NGX, NGY
      CHARACTER*1       SCHEME
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,7), RHS(LDA)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDY, PDEF
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS
      INTEGER           IERROR, MINLDA, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          D03EEZ
C     .. Executable Statements ..
C
C     Check for invalid parameters first --
C
      EPS = X02AJF()
      IERROR = 0
      MINLDA = NGX*NGY
      IF (XMIN.GE.XMAX) THEN
         IERROR = 1
         WRITE (REC,FMT=99999) XMIN, XMAX
         NREC = 2
      ELSE IF (YMIN.GE.YMAX) THEN
         IERROR = 1
         WRITE (REC,FMT=99998) YMIN, YMAX
         NREC = 2
      ELSE IF (NGX.LT.3) THEN
         IERROR = 1
         WRITE (REC,FMT=99996) NGX
         NREC = 1
      ELSE IF (NGY.LT.3) THEN
         IERROR = 1
         WRITE (REC,FMT=99995) NGY
         NREC = 1
      ELSE IF (LDA.LT.MINLDA) THEN
         IERROR = 1
         WRITE (REC,FMT=99997) MINLDA, LDA
         NREC = 2
      ELSE IF (SCHEME.NE.'C' .AND. SCHEME.NE.'c' .AND. SCHEME.NE.
     *         'U' .AND. SCHEME.NE.'u') THEN
         IERROR = 1
         WRITE (REC,FMT=99994) SCHEME
         NREC = 2
      END IF
      IF (IERROR.EQ.0) THEN
         CALL D03EEZ(NGX,NGY,LDA,A,RHS,XMIN,XMAX,YMIN,YMAX,PDEF,BNDY,
     *               SCHEME,IERROR,IFAIL)
         NREC = 0
         IF (IERROR.EQ.4) THEN
            REC(1) = ' ** Equation not elliptic at some point'
            NREC = 1
         ELSE IF (IERROR.EQ.5) THEN
            REC(1) = ' ** There is no unique solution with Neumann'//
     *               ' Boundary conditions'
            NREC = 1
         ELSE IF (IERROR.EQ.6) THEN
            REC(1) =
     *          ' ** The linear equations were not diagonally dominated'
            NREC = 1
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** On entry, XMAX must be greater than XMIN:',/' ** XM',
     *       'IN = ',1P,D13.5,',  XMAX = ',D13.5)
99998 FORMAT (' ** On entry, YMAX must be greater than YMIN:',/' ** YM',
     *       'IN = ',1P,D13.5,',  YMAX = ',D13.5)
99997 FORMAT (' ** On entry LDA must be at least NGX*NGY:',/' ** NGX*N',
     *       'GY = ',I16,',  LDA = ',I16)
99996 FORMAT (' ** On entry, NGX is less than 3: NGX = ',I16)
99995 FORMAT (' ** On entry, NGY is less than 3: NGY = ',I16)
99994 FORMAT (' ** On entry, SCHEME should be one of ''C'', ''c'',''U',
     *       ''' or ''u''',/' ** SCHEME = ''',A1,'''')
      END
