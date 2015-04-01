      DOUBLE PRECISION FUNCTION S21BCF(X,Y,Z,IFAIL)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C
C     ******************************************************************
C
C      Calculates an approximate value for the elliptic integral of
C      the second kind, Rd(X,Y,Z), as defined by B.C.Carlson.
C       - Special Functions of Applied Maths - Academic Press(1977).
C      The algorithm is also due to Carlson.
C
C      Rd(X,Y,Z)=integral(zero to infinity) of
C          1.5/sqrt((t+X)*(t+Y)*(t+Z)**3) dt
C
C      for X.ge.0.0 , Y.ge.0.0, Z.gt.0.0 and at most one of X
C      and Y equal to 0.0.
C
C     ******************************************************************
C     Precision-dependent constant
C
C        ACC  Controls final truncation accuracy. An accuracy
C                 close to full machine precision can be obtained
C                 for all legal arguments if ACC is chosen so that
C           3.0*ACC**6/(1.0-ACC)**1.5.le.X02AJF/X02BHF
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
      PARAMETER                        (SRNAME='S21BCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X, Y, Z
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 ACC, C1, C2, C3, C4, C5, C6, CXN,
     *                                 CYN, CZN, EA, EB, EC, ED, EF,
     *                                 LAMDA, LOLIM, MU, POW4, RD, RTX,
     *                                 RTY, RTZ, SIGMA, UPLIM, XN, YN,
     *                                 ZN
      INTEGER                          IND, NREC
C     .. Local Arrays ..
      CHARACTER*80                     P01REC(3)
C     .. External Functions ..
      DOUBLE PRECISION                 X02AMF
      INTEGER                          P01ABF
      EXTERNAL                         X02AMF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MAX, MIN, SQRT
C     .. Data statements ..
C      * EXPANSION (DATA) *
C08   DATA                             ACC/3.8D-2/
C09   DATA                             ACC/2.6D-2/
C12   DATA                             ACC/8.3D-3/
C15   DATA                             ACC/2.6D-3/
      DATA                             ACC/1.2D-3/
C19   DATA                             ACC/5.6D-4/
C     .. Executable Statements ..
C
      IND = 0
      NREC = 0
      RD = 0.0D0
      LOLIM = 2.0D0*X02AMF()**0.66667D0
      UPLIM = (0.1D0*ACC/X02AMF())**0.66667D0
C  Test for valid arguments
      IF (X.LT.0.0D0 .OR. Y.LT.0.0D0 .OR. X+Y.EQ.0.0D0) THEN
         IND = 1
         NREC = 1
         WRITE (P01REC,FMT=99999) X, Y
      ELSE IF (Z.LE.0.0D0) THEN
         IND = 2
         NREC = 1
         WRITE (P01REC,FMT=99998) Z
      ELSE IF (MIN(X+Y,Z).LT.LOLIM) THEN
         IND = 3
         NREC = 3
         WRITE (P01REC,FMT=99997) LOLIM, X, Y, Z
      ELSE IF (MAX(X,Y,Z).GE.UPLIM) THEN
         IND = 4
         NREC = 2
         WRITE (P01REC,FMT=99996) X, Y, Z
      ELSE
C
C     Valid call
         XN = X
         YN = Y
         ZN = Z
         SIGMA = 0.0D0
         POW4 = 1.0D0
C
C     Main recursion
   20    MU = (XN+YN+3.0D0*ZN)*0.2D0
         CXN = (MU-XN)/MU
         CYN = (MU-YN)/MU
         CZN = (MU-ZN)/MU
         IF (MAX(ABS(CXN),ABS(CYN),ABS(CZN)).GE.ACC) THEN
            RTX = SQRT(XN)
            RTY = SQRT(YN)
            RTZ = SQRT(ZN)
            LAMDA = RTX*RTY + RTY*RTZ + RTZ*RTX
            ZN = ZN + LAMDA
            SIGMA = SIGMA + POW4/(RTZ*ZN)
            POW4 = POW4*0.25D0
            ZN = ZN*0.25D0
            XN = (XN+LAMDA)*0.25D0
            YN = (YN+LAMDA)*0.25D0
            GO TO 20
         END IF
C
C     Power series evaluation
         C1 = 3.0D0/14.0D0
         C2 = 9.0D0/88.0D0
         C3 = 9.0D0/52.0D0
         C4 = 1.0D0/6.0D0
         C5 = 9.0D0/22.0D0
         C6 = 3.0D0/26.0D0
         EA = CXN*CYN
         EB = CZN*CZN
         EC = 3.0D0*EA - 8.0D0*EB
         ED = EA - EB
         EF = EA - 6.0D0*EB
         RD = 3.0D0*SIGMA + POW4*(1.0D0+EF*(C2*EF-C1-C3*EC*CZN)
     *        +CZN*(C4*EC+CZN*(CZN*C6*EA-C5*ED)))/(MU*SQRT(MU))
      END IF
      S21BCF = RD
      IFAIL = P01ABF(IFAIL,IND,SRNAME,0,P01REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, X or Y is less than zero:',1P,'X =',
     *       D13.5,', Y =',D13.5)
99998 FORMAT (1X,'** On entry, Z is .le. zero:',1P,'Z =',D13.5)
99997 FORMAT (1X,'** On entry, either (X+Y) or Z is less than ',1P,
     *       D13.5,':',/4X,'X =',D13.5,', Y =',D13.5,', Z =',D13.5,/4X,
     *       'There is a danger of overflow')
99996 FORMAT (1X,'** On entry, at least one of X, Y and Z is greater t',
     *       'han ',1P,D13.5,':',/4X,'X =',D13.5,', Y =',D13.5,', Z =',
     *       D13.5,/4X,'There is a danger of underflow')
      END
