      SUBROUTINE C02AHF(AR,AI,BR,BI,CR,CI,ZSM,ZLG,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     C02AHF DETERMINES THE ROOTS OF THE QUADRATIC EQUATION
C        A*Z**2 + B*Z + C = 0
C     WHERE A = CMPLX(AR,AI), B = CMPLX(BR,BI) AND C = CMPLX(CR,CI)
C     ARE COMPLEX COEFFICIENTS, AND ZSM AND ZLG ARE THE SMALLEST AND
C     LARGEST ROOT IN MAGNITUDE RESPECTIVELY.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C02AHF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  AI, AR, BI, BR, CI, CR
      INTEGER           IFAIL
C     .. Array Arguments ..
      DOUBLE PRECISION  ZLG(2), ZSM(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  APRIME, BPRIME, CPRIME, FINITY
      INTEGER           EMAXM1, EMINM1, EXPDEP, IER, NREC
      LOGICAL           OVFLOW
C     .. Local Arrays ..
      CHARACTER*80      REC(6)
C     .. External Functions ..
      DOUBLE PRECISION  X02ALF
      INTEGER           P01ABF, X02BJF, X02BKF, X02BLF
      EXTERNAL          X02ALF, P01ABF, X02BJF, X02BKF, X02BLF
C     .. External Subroutines ..
      EXTERNAL          C02AFW, C02AHZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
      NREC = 2
      FINITY = X02ALF()
      IER = 0
      IF ((AR.EQ.ZERO) .AND. (AI.EQ.ZERO)) THEN
         IER = 1
         WRITE (REC,FMT=99999)
         NREC = 1
         IF ((BR.EQ.ZERO) .AND. (BI.EQ.ZERO)) THEN
            IER = 2
            WRITE (REC,FMT=99998)
            NREC = 1
            ZLG(1) = FINITY
            ZLG(2) = ZERO
            ZSM(1) = FINITY
            ZSM(2) = ZERO
         ELSE
            ZLG(1) = FINITY
            ZLG(2) = ZERO
            CALL C02AFW(-CR,-CI,BR,BI,ZSM(1),ZSM(2),OVFLOW)
            IF (OVFLOW) THEN
               IER = 3
               NREC = 3
               WRITE (REC,FMT=99997) AR, CR, BR, AI, CI, BI
               ZSM(1) = FINITY
               ZSM(2) = ZERO
            END IF
         END IF
      ELSE IF ((CR.EQ.ZERO) .AND. (CI.EQ.ZERO)) THEN
         ZSM(1) = ZERO
         ZSM(2) = ZERO
         CALL C02AFW(-BR,-BI,AR,AI,ZLG(1),ZLG(2),OVFLOW)
         IF (OVFLOW) THEN
            IER = 4
            NREC = 3
            WRITE (REC,FMT=99996) CR, BR, AR, CI, BI, AI
            ZLG(1) = FINITY
            ZLG(2) = ZERO
         END IF
      ELSE
         EXPDEP = X02BJF() + 1
         EMINM1 = X02BKF() - 1
         EMAXM1 = X02BLF() - 1
         CALL C02AHZ(EXPDEP,EMINM1,EMAXM1,FINITY,AR,AI,BR,BI,CR,CI,ZSM,
     *               ZLG,IER)
         IF (IER.EQ.0) THEN
            IFAIL = 0
            GO TO 20
         ELSE
            NREC = 6
            APRIME = MAX(ABS(AR),ABS(AI))
            BPRIME = MAX(ABS(BR),ABS(BI))
            CPRIME = MAX(ABS(CR),ABS(CI))
            WRITE (REC,FMT=99995) BPRIME, APRIME, CPRIME, BR, BI, AR, AI
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
   20 RETURN
C
99999 FORMAT (' ** On entry, (AR,AI).eq.(0,0)')
99998 FORMAT (' ** On entry, (AR,AI).eq.(0,0) and (BR,BI).eq.(0,0)')
99997 FORMAT (' ** On entry, (AR,AI).eq.(0,0) and the root -(CR,CI)/(B',
     *       'R,BI) overflows:',/'    AR = ',1P,D13.5,' CR = ',1P,D13.5,
     *       ' BR = ',1P,D13.5,/'    AI = ',1P,D13.5,' CI = ',1P,D13.5,
     *       ' BI = ',1P,D13.5)
99996 FORMAT (' ** On entry, (CR,CI).eq.(0,0) and the root -(BR,BI)/(A',
     *       'R,AI) overflows:',/'    CR = ',1P,D13.5,' BR = ',1P,D13.5,
     *       ' AR = ',1P,D13.5,/'    CI = ',1P,D13.5,' BI = ',1P,D13.5,
     *       ' AI = ',1P,D13.5)
99995 FORMAT (' ** On entry, B'' is so large that B''**2 is indistingu',
     *       'ishable',/'    from (B''**2 - 4*A''*C'') and the root -(',
     *       'BR,BI)/(AR,AI) overflows:',/'    B'' = MAX(ABS(BR),ABS(B',
     *       'I)) = ',1P,D13.5,/'    A'' = MAX(ABS(AR),ABS(AI)) = ',1P,
     *       D13.5,/'    C'' = MAX(ABS(CR),ABS(CI)) = ',1P,D13.5,/'   ',
     *       ' BR = ',1P,D13.5,' BI = ',1P,D13.5,' AR = ',1P,D13.5,' A',
     *       'I = ',1P,D13.5)
      END
