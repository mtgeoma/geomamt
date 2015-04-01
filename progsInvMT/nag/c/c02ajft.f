      SUBROUTINE C02AJF(A,B,C,ZSM,ZLG,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     C02AJF DETERMINES THE ROOTS OF THE QUADRATIC EQUATION
C        A*Z**2 + B*Z + C = 0
C     WHERE A, B AND C ARE REAL COEFFICIENTS, AND ZSM AND ZLG
C     ARE THE SMALLEST AND LARGEST ROOT IN MAGNITUDE RESPECTIVELY.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C02AJF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, C
      INTEGER           IFAIL
C     .. Array Arguments ..
      DOUBLE PRECISION  ZLG(2), ZSM(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  FINITY
      INTEGER           EMAXM1, EMINM1, EXPDEP, IER, NREC
      LOGICAL           OVFLOW
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  F06BLF, X02ALF
      INTEGER           P01ABF, X02BJF, X02BKF, X02BLF
      EXTERNAL          F06BLF, X02ALF, P01ABF, X02BJF, X02BKF, X02BLF
C     .. External Subroutines ..
      EXTERNAL          C02AJZ
C     .. Executable Statements ..
      NREC = 2
      FINITY = X02ALF()
      IER = 0
      IF (A.EQ.ZERO) THEN
         IER = 1
         WRITE (REC,FMT=99999)
         NREC = 1
         IF (B.EQ.ZERO) THEN
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
            ZSM(1) = F06BLF(-C,B,OVFLOW)
            ZSM(2) = ZERO
            IF (OVFLOW) THEN
               IER = 3
               WRITE (REC,FMT=99997) A, C, B
               NREC = 2
               ZSM(1) = FINITY
            END IF
         END IF
      ELSE IF (C.EQ.ZERO) THEN
         ZSM(1) = ZERO
         ZSM(2) = ZERO
         ZLG(1) = F06BLF(-B,A,OVFLOW)
         ZLG(2) = ZERO
         IF (OVFLOW) THEN
            IER = 4
            WRITE (REC,FMT=99996) C, B, A
            ZLG(1) = FINITY
         END IF
      ELSE
         EXPDEP = X02BJF() + 1
         EMINM1 = X02BKF() - 1
         EMAXM1 = X02BLF() - 1
         CALL C02AJZ(EXPDEP,EMINM1,EMAXM1,FINITY,A,B,C,ZSM,ZLG,IER)
         IF (IER.EQ.0) THEN
            IFAIL = 0
            GO TO 20
         ELSE
            NREC = 3
            WRITE (REC,FMT=99995) B, A, C
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
   20 RETURN
C
99999 FORMAT (' ** On entry, A.eq.0')
99998 FORMAT (' ** On entry, A.eq.0 and B.eq.0')
99997 FORMAT (' ** On entry, A.eq.0 and the root -C/B overflows:',/'  ',
     *       '  A = ',1P,D13.5,' C = ',1P,D13.5,' B = ',1P,D13.5)
99996 FORMAT (' ** On entry, C.eq.0 and the root -B/A overflows:',/'  ',
     *       '  C = ',1P,D13.5,' B = ',1P,D13.5,' A = ',1P,D13.5)
99995 FORMAT (' ** On entry, B is so large that B**2 is indistinguisha',
     *       'ble',/'    from (B**2 - 4*A*C) and the root -B/A overflo',
     *       'ws:',/'    B = ',1P,D13.5,' A = ',1P,D13.5,' C = ',1P,
     *       D13.5)
      END
