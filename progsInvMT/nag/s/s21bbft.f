      DOUBLE PRECISION FUNCTION S21BBF(X,Y,Z,IFAIL)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C
C     ******************************************************************
C
C      Calculates an approximate value for the elliptic integral of
C      the first kind, Rf(X,Y,Z), as defined by B.C.Carlson.
C      - Special Functions of Applied Maths - Academic Press(1977).
C      The algorithm is also due to Carlson but has been recoded by
C      J.L.Schonfelder to avoid intermediate under and overflows.
C
C      Rf(X,Y,Z)=integral(zero to infinity) of
C          0.5/sqrt((t+X)*(t+Y)*(t+Z)) dt
C
C      for X,Y,Z all .ge. zero and at most one equal to zero.
C
C     ******************************************************************
C     Precision-dependent constant
C
C        ACC  Controls final truncation accuracy. An accuracy
C                 close to full machine precision can be obtained
C                 for all legal arguments if ACC is chosen so that
C           0.25*ACC**6/(1.0-ACC).le.X02AJF/X02BHF
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
      PARAMETER                        (SRNAME='S21BBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X, Y, Z
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 ACC, C1, C2, C3, CXN, CYN, CZN,
     *                                 E2, E3, LAMDA, LOLIM, MU, RF,
     *                                 RSCALE, RTX, RTY, RTZ, UPLIM, XN,
     *                                 YN, ZN
      INTEGER                          IND, NREC
C     .. Local Arrays ..
      CHARACTER*80                     P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION                 X02AKF, X02ALF
      INTEGER                          P01ABF
      EXTERNAL                         X02AKF, X02ALF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        MAX, SQRT
C     .. Data statements ..
C      * EXPANSION (DATA) *
C08   DATA                             ACC/5.7D-2/
C09   DATA                             ACC/3.9D-2/
C12   DATA                             ACC/1.2D-2/
C15   DATA                             ACC/3.9D-3/
      DATA                             ACC/1.8D-3/
C19   DATA                             ACC/8.5D-4/
C     .. Executable Statements ..
C
      IND = 0
      NREC = 0
      RF = 0.0D0
C  Order X,Y,Z into XN,YN,ZN such that XN.le.YN.le.ZN
      IF (X.GT.Y) THEN
         XN = Y
         YN = X
      ELSE
         XN = X
         YN = Y
      END IF
      ZN = Z
      IF (YN.GT.ZN) THEN
         ZN = YN
         YN = Z
         IF (XN.GT.YN) THEN
            YN = XN
            XN = Z
         END IF
      END IF
C
C  Test for valid arguments
      IF (XN.LT.0.0D0) THEN
         IND = 1
         NREC = 2
         WRITE (P01REC,FMT=99999) X, Y, Z
      ELSE IF (YN.EQ.0.0D0) THEN
         IND = 2
         NREC = 2
         WRITE (P01REC,FMT=99998) X, Y, Z
      ELSE
C
C     Valid call
         RSCALE = 1.0D0
         LOLIM = 16.0D0*X02AKF()
         UPLIM = 0.0625D0*X02ALF()
C
C     For extreme arguments scale to avoid under and overflows
         IF (ZN.GT.UPLIM) THEN
            RSCALE = 0.25D0
            ZN = ZN*0.0625D0
            IF (YN.GT.LOLIM) THEN
               YN = YN*0.0625D0
               IF (XN.GT.LOLIM) THEN
                  XN = XN*0.0625D0
               ELSE
                  RTZ = SQRT(ZN)
                  RTY = SQRT(YN)
                  LAMDA = RTZ*RTY + ((RTZ+RTY)*0.25D0)*SQRT(XN)
                  XN = LAMDA*0.25D0
                  YN = (YN+LAMDA)*0.25D0
                  ZN = (ZN+LAMDA)*0.25D0
               END IF
            ELSE
               LAMDA = (SQRT(XN)+SQRT(YN))*(SQRT(ZN)*0.25D0)
               XN = LAMDA*0.25D0
               YN = XN
               ZN = (ZN+LAMDA)*0.25D0
            END IF
         ELSE IF (ZN.LE.LOLIM) THEN
            RSCALE = 4.0D0
            XN = XN*16.0D0
            YN = YN*16.0D0
            ZN = ZN*16.0D0
         END IF
C
C     Main recursion
   20    MU = (XN+YN+ZN)/3.0D0
         CZN = 2.0D0 - (ZN+MU)/MU
         CXN = 2.0D0 - (XN+MU)/MU
         IF (MAX(CXN,-CZN).GT.ACC) THEN
            RTX = SQRT(XN)
            RTY = SQRT(YN)
            RTZ = SQRT(ZN)
            LAMDA = RTZ*(RTX+RTY) + RTX*RTY
            XN = (XN+LAMDA)*0.25D0
            YN = (YN+LAMDA)*0.25D0
            ZN = (ZN+LAMDA)*0.25D0
            GO TO 20
         END IF
C
C     Power series expansion
         C1 = 1.0D0/24.0D0
         C2 = 3.0D0/44.0D0
         C3 = 1.0D0/14.0D0
         CYN = -CXN - CZN
         E2 = CXN*CYN - CZN*CZN
         E3 = CXN*CZN*CYN
         RF = RSCALE*(1.0D0+(C1*E2-0.1D0-C2*E3)*E2+C3*E3)/SQRT(MU)
      END IF
      S21BBF = RF
      IFAIL = P01ABF(IFAIL,IND,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, at least one of X, Y and Z is negative:',
     *       /4X,'X =',1P,D13.5,', Y =',D13.5,', Z =',D13.5)
99998 FORMAT (1X,'** On entry, two or more of X, Y and Z are zero:',/4X,
     *       'X =',1P,D13.5,', Y =',D13.5,', Z =',D13.5)
      END
