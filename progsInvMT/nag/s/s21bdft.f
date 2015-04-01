      DOUBLE PRECISION FUNCTION S21BDF(X,Y,Z,R,IFAIL)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C
C     ******************************************************************
C
C      Calculates an approximate value for the elliptic integral of
C      the third kind, Rj(X,Y,Z,R), as defined by B.C.Carlson.
C      - Special Functions of Applied Maths - Academic Press(1977).
C      The algorithm is also due to Carlson
C
C      Rj(X,Y,Z,R)=integral(zero to infinity) of
C          1.5/(sqrt((t+X)*(t+Y)*(t+Z))*(t+R)) dt
C
C      for X.ge.0.0 , Y.ge.0.0 , Z.ge.0.0
C         X+Y.gt.0.0
C         Y+Z.gt.0.0
C         Z+X.gt.0.0
C         R.ne.0.0
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
C     activate the statement contained in a comment beginning  CDD,
C     where  DD  is the approximate number of significant decimal
C     digits represented by the machine.
C     * EXPANSION (DATA) *
C     ******************************************************************
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='S21BDF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 R, X, Y, Z
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 A, ACC, AN, B, BN, C1, C2, C3,
     *                                 C4, C5, C6, CRN, CXN, CYN, CZN,
     *                                 E2, E3, EA, EB, EC, LAMDA, LOLIM,
     *                                 MU, POW4, RCX, RHO, RJ, RN, RTX,
     *                                 RTY, RTZ, S1, S2, S3, SIGMA, TAU,
     *                                 UPLIM, XN, YN, YY, ZN
      INTEGER                          IND, JFAIL, NREC
C     .. Local Arrays ..
      CHARACTER*80                     P01REC(3)
C     .. External Functions ..
      DOUBLE PRECISION                 S21BAF, S21BBF, X02AMF
      INTEGER                          P01ABF
      EXTERNAL                         S21BAF, S21BBF, X02AMF, P01ABF
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
      RJ = 0.0D0
      LOLIM = X02AMF()**0.33333D0
      UPLIM = (1.0D0/X02AMF())**0.33333D0*0.396D0
      IND = 0
      NREC = 0
C  Test for valid arguments
      IF (X.LT.0.0D0 .OR. Y.LT.0.0D0 .OR. Z.LT.0.0D0) THEN
         IND = 1
         NREC = 2
         WRITE (P01REC,FMT=99999) X, Y, Z
      ELSE IF (X+Y.EQ.0.0D0 .OR. Y+Z.EQ.0.0D0 .OR. Z+X.EQ.0.0D0) THEN
         IND = 1
         NREC = 2
         WRITE (P01REC,FMT=99998) X, Y, Z
      ELSE IF (R.EQ.0.0D0) THEN
         IND = 2
         NREC = 1
         WRITE (P01REC,FMT=99997)
      ELSE IF (MIN(X+Y,Y+Z,Z+X,ABS(R)).LT.LOLIM) THEN
         IND = 3
         NREC = 3
         WRITE (P01REC,FMT=99996) LOLIM, R, X, Y, Z
      ELSE IF (MAX(X,Y,Z,ABS(R)).GT.UPLIM) THEN
         IND = 4
         NREC = 2
         WRITE (P01REC,FMT=99995) UPLIM, R, X, Y, Z
      ELSE
C     Valid call
         IF (R.GT.0.0D0) THEN
            XN = X
            YN = Y
            ZN = Z
            RN = R
         ELSE
C        Order X, Y and Z, and transform to positive R
            XN = MIN(X,Y)
            YY = MAX(X,Y)
            ZN = MAX(YY,Z)
            YY = MIN(YY,Z)
            YN = MAX(XN,YY)
            XN = MIN(XN,YY)
            A = 1.0D0/(YN-R)
            B = (ZN-YN)*A*(YN-XN)
            RN = YN + B
            RHO = XN*ZN/YN
            TAU = R*RN/YN
            JFAIL = 0
            RCX = S21BAF(RHO,TAU,JFAIL)
         END IF
         SIGMA = 0.0D0
         POW4 = 1.0D0
C
C     Main recursion
   20    MU = (XN+YN+ZN+2.0D0*RN)*0.2D0
         CXN = (MU-XN)/MU
         CYN = (MU-YN)/MU
         CZN = (MU-ZN)/MU
         CRN = (MU-RN)/MU
         IF (MAX(ABS(CXN),ABS(CYN),ABS(CZN),ABS(CRN)).GE.ACC) THEN
            RTX = SQRT(XN)
            RTY = SQRT(YN)
            RTZ = SQRT(ZN)
            LAMDA = RTX*RTY + RTY*RTZ + RTZ*RTX
            AN = RN*(RTX+RTY+RTZ) + RTX*RTY*RTZ
            AN = AN*AN
            BN = RN*(RN+LAMDA)**2
            JFAIL = 0
            SIGMA = SIGMA + POW4*S21BAF(AN,BN,JFAIL)
            POW4 = POW4*0.25D0
            XN = (XN+LAMDA)*0.25D0
            YN = (YN+LAMDA)*0.25D0
            ZN = (ZN+LAMDA)*0.25D0
            RN = (RN+LAMDA)*0.25D0
            GO TO 20
         END IF
C
C     Power series evaluation
   40    C1 = 3.0D0/14.0D0
         C2 = 9.0D0/88.0D0
         C3 = 9.0D0/52.0D0
         C4 = 1.0D0/3.0D0
         C5 = 3.0D0/22.0D0
         C6 = 3.0D0/26.0D0
         EA = CXN*CYN + CYN*CZN + CZN*CXN
         EB = CXN*CYN*CZN
         EC = CRN*CRN
         E2 = EA - 3.0D0*EC
         E3 = EB + 2.0D0*CRN*(EA-EC)
         S1 = 1.0D0 + E2*(C2*E2-C1-C3*E3)
         S2 = EB*(0.5D0*C4+CRN*(CRN*C6-C5*2.0D0))
         S3 = CRN*EA*(C4-CRN*C5) - C4*CRN*EC
         RJ = 3.0D0*SIGMA + POW4*(S1+S2+S3)/(MU*SQRT(MU))
         IF (R.LT.0.0D0) RJ = A*(B*RJ+3.0D0*(RCX-S21BBF(XN,YN,ZN,JFAIL))
     *                        )
      END IF
C
      S21BDF = RJ
      IFAIL = P01ABF(IFAIL,IND,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, one of X, Y or Z is less than zero:',/4X,
     *       1P,'X =',D13.5,', Y =',D13.5,', Z =',D13.5)
99998 FORMAT (1X,'** On entry, more than one of X, Y and Z are equal t',
     *       'o zero:',/4X,1P,'X =',D13.5,', Y =',D13.5,', Z =',D13.5)
99997 FORMAT (1X,'** On entry, R is equal to zero')
99996 FORMAT (1X,'** On entry, either abs(R), or the sum of any two of',
     *       ' X, Y and Z, is less',/4X,'than ',1P,D13.5,'. There is a',
     *       ' danger of overflow.',/4X,1P,'R =',D13.5,', X =',D13.5,
     *       ', Y =',D13.5,', Z =',D13.5)
99995 FORMAT (1X,'** On entry, at least one of abs(R), X, Y and Z is g',
     *       'reater than',1P,D13.5,' :',/4X,1P,'R =',D13.5,', X =',
     *       D13.5,', Y =',D13.5,', Z =',D13.5)
      END
