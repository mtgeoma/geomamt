      DOUBLE PRECISION FUNCTION S21BAF(X,Y,IFAIL)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C
C     ******************************************************************
C
C     Calculates an approximate value for the function
C     Rc(x,y) related to the elliptic integrals.
C     The function is as defined by  Carlson - Special
C           Functions of Applied Mathematics - Academic Press(1977)
C
C     Rc(x,y)=integral(zero to infinity) of
C               0.5/(sqrt(t+x)*(t+y)) dt
C
C     subject to restrictions
C          x.ge.0.0 and y.ne.0.0
C     ******************************************************************
C     Precision-dependent constant
C
C        ACC  Controls final truncation accuracy. An accuracy
C                 close to full machine precision can be obtained
C                 for all legal arguments if ACC is chosen so that
C           16.0*ACC**6/(1.0-2.0*ACC).le.X02AJF/X02BHF
C
C     To select the correct value for a particular machine-range,
C     activate the statement contained in a comment beginning  CDD ,
C     where  DD  is the approximate number of significant decimal
C     digits represented by the machine.
C     * EXPANSION (DATA) *
C     ******************************************************************
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='S21BAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X, Y
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 ACC, C3, C5, LAMDA, LOLIM, MU,
     *                                 RC, RSCALE, SN, TMP, UPLIM, XN,
     *                                 YN
      INTEGER                          IND, NREC
C     .. Local Arrays ..
      CHARACTER*80                     P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 X02AKF, X02ALF
      INTEGER                          P01ABF
      EXTERNAL                         X02AKF, X02ALF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, SQRT
C     .. Data statements ..
C      * EXPANSION (DATA) *
C08   DATA                             ACC/2.8D-2/
C09   DATA                             ACC/1.9D-2/
C12   DATA                             ACC/6.2D-3/
C15   DATA                             ACC/1.9D-3/
      DATA                             ACC/9.2D-4/
C19   DATA                             ACC/4.2D-4/
C     .. Executable Statements ..
      RC = 0.0D0
      IND = 0
      NREC = 1
C  Test for argument validity
      IF (X.LT.0.0D0) THEN
         IND = 1
         NREC = 1
         WRITE (P01REC,FMT=99999) X
      ELSE IF (Y.EQ.0.0D0) THEN
         IND = 2
         NREC = 1
         WRITE (P01REC,FMT=99998)
      ELSE
C
C     Valid call
         IF (Y.GT.0.0D0) THEN
            XN = X
            YN = Y
            RSCALE = 1.0D0
         ELSE
            XN = X - Y
            YN = -Y
            RSCALE = SQRT(X)/SQRT(XN)
         END IF
         LOLIM = 16.0D0*X02AKF()
         UPLIM = 0.0625D0*X02ALF()
C
C     For arguments in extreme regions apply
C     scaling to avoid intermediate underflows
C     and overflows
         IF (XN.LE.LOLIM .AND. YN.LE.LOLIM) THEN
            RSCALE = 4.0D0*RSCALE
            XN = XN*16.0D0
            YN = YN*16.0D0
         ELSE IF (XN.GE.UPLIM .OR. YN.GE.UPLIM) THEN
            RSCALE = 0.25D0*RSCALE
            IF (XN.LT.LOLIM) THEN
               XN = YN*0.015625D0
               YN = 2.0D0*XN
            ELSE IF (YN.LT.LOLIM) THEN
               TMP = XN*0.015625D0
               YN = 0.03125D0*SQRT(XN*YN)
               XN = TMP
            ELSE
               XN = XN*0.0625D0
               YN = YN*0.0625D0
            END IF
         END IF
C
C     Main recursion
   20    MU = (XN+2.0D0*YN)/3.0D0
         SN = (YN+MU)/MU - 2.0D0
         IF (ABS(SN).LT.ACC) GO TO 40
         LAMDA = 2.0D0*SQRT(XN)*SQRT(YN) + YN
         XN = (XN+LAMDA)*0.25D0
         YN = (YN+LAMDA)*0.25D0
         GO TO 20
C
C     Taylor expansion
   40    C3 = 1.0D0/7.0D0
         C5 = 9.0D0/22.0D0
         RC = RSCALE*(1.0D0+SN*SN*(0.3D0+SN*(C3+SN*(0.375D0+SN*C5))))
     *        /SQRT(MU)
      END IF
C
      S21BAF = RC
      IFAIL = P01ABF(IFAIL,IND,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, X is less than zero: X =',1P,D13.5)
99998 FORMAT (1X,'** On entry, Y is equal to zero')
      END
