      DOUBLE PRECISION FUNCTION G01EMF(Q,V,IR,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     COMPUTES THE PROBABILITY INTEGRAL FROM 'ZERO' TO 'Q'.
C     OF THE STUDENTIZED RANGE DISTRIBUTION
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      INTEGER                          INF, LIW, LW
      PARAMETER                        (SRNAME='G01EMF',INF=2,LIW=250,
     *                                 LW=1000)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 Q, V
      INTEGER                          IFAIL, IR
C     .. Scalars in Common ..
      DOUBLE PRECISION                 CC, CEPS, CQ, CR, CV, UFLOW
      INTEGER                          ICR
C     .. Local Scalars ..
      DOUBLE PRECISION                 ABSACC, ABSERR, BOUND, EPS,
     *                                 EPSREL, LIM, PEE, R1, RESULT, V2,
     *                                 YA, YB
      INTEGER                          IERROR, IFAULT, NPTS, NREC
C     .. Local Arrays ..
      DOUBLE PRECISION                 W(LW)
      INTEGER                          IW(LIW)
      CHARACTER*80                     P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION                 G01EMW, G01EMX, G01EMY, G01EMZ,
     *                                 S14ABF, S15ABF, X02AMF
      INTEGER                          P01ABF
      EXTERNAL                         G01EMW, G01EMX, G01EMY, G01EMZ,
     *                                 S14ABF, S15ABF, X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL                         D01AMF, D01DAF
C     .. Intrinsic Functions ..
      INTRINSIC                        EXP, LOG, MAX, MIN, DBLE
C     .. Common blocks ..
      COMMON                           /AG01EM/CQ, CR, CV, CC, CEPS,
     *                                 UFLOW, ICR
C     .. Executable Statements ..
C
      RESULT = 0.0D0
      NREC = 1
      IERROR = 1
      IF (V.LT.1.0D0) THEN
         WRITE (P01REC,FMT=99999) V
      ELSE IF (IR.LT.2) THEN
         WRITE (P01REC,FMT=99998) IR
      ELSE IF (Q.LE.0.0D0) THEN
         WRITE (P01REC,FMT=99997) Q
      ELSE
C
C        Set up bounds and call integrator
C
         CQ = Q
         IERROR = 0
         ICR = IR - 1
         CR = DBLE(IR)
         ABSACC = 0.000005D0
         IF (V.LE.2000) THEN
            IF (V.LE.5.0D0) THEN
               YA = 0.0D0
               YB = 2.61D0 + EXP(1.64D0-0.485D0*V)
            ELSE IF (V.LE.25.0D0) THEN
               YA = 0.74D0 - EXP(-0.14D0-0.0282D0*V)
               YB = 1.6D0 + EXP(0.7D0-0.0778D0*V)
            ELSE IF (V.LE.200.0D0) THEN
               YA = 0.74D0 - EXP(-0.5D0-0.0147D0*V)
               YB = 1.25D0 + EXP(-0.15D0-0.013D0*V)
            ELSE
               YA = 0.94D0 - EXP(-1.3D0-0.000922D0*V)
               YB = 1.08D0 + EXP(-1.2D0-0.0015D0*V)
            END IF
            CV = V
            R1 = 1.0D0/(DBLE(ICR))
            EPS = X02AMF()
            UFLOW = LOG(EPS)
            LIM = 1.5D0/Q*(EPS**R1)
            YA = MAX(YA,LIM)
            V2 = V/2.0D0
            IFAULT = 0
            CC = (V2)*LOG(V) - S14ABF(V2,IFAULT) - (V2-1.0D0)*LOG(2.0D0)
            PEE = 2.0D0*S15ABF(0.5D0*Q*YB,IFAULT) - 1.0D0
            EPS = 1.0D-07*PEE
            CEPS = EPS**R1
            IFAULT = 1
            CALL D01DAF(YA,YB,G01EMX,G01EMW,G01EMY,ABSACC,RESULT,NPTS,
     *                  IFAULT)
         ELSE
C
C           Use assymptotic integral
C
            CV = 0.0D0
            CC = 0.0D0
            CEPS = 0.0D0
            UFLOW = 0.0D0
            EPSREL = 0.0D0
            IFAULT = 1
            CALL D01AMF(G01EMZ,BOUND,INF,ABSACC,EPSREL,RESULT,ABSERR,W,
     *                  LW,IW,LIW,IFAULT)
         END IF
         IF (IFAULT.NE.0) THEN
            NREC = 2
            IERROR = 2
            WRITE (P01REC,FMT=99996)
         END IF
      END IF
      G01EMF = MIN(1.0D0,MAX(0.0D0,RESULT))
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, V.lt.1.0 : V = ',D13.5)
99998 FORMAT (' ** On entry, IR.lt.2 : IR = ',I16)
99997 FORMAT (' ** On entry, Q.le.0.0 : Q = ',D13.5)
99996 FORMAT (' ** Warning - There is some doubt as to whether',/'    ',
     *       '          full accuracy has been achieved. ')
      END
