      COMPLEX*16  FUNCTION S01EAF(Z,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     Returns exp(Z) for complex Z.
C     .. Parameters ..
      DOUBLE PRECISION            ONE, ZERO
      PARAMETER                   (ONE=1.0D0,ZERO=0.0D0)
      CHARACTER*6                 SRNAME
      PARAMETER                   (SRNAME='S01EAF')
C     .. Scalar Arguments ..
      COMPLEX*16                  Z
      INTEGER                     IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION            COSY, EXPX, LNSAFE, RECEPS, RESI,
     *                            RESR, RTSAFS, SAFE, SAFSIN, SINY, X,
     *                            XPLNCY, XPLNSY, Y
      INTEGER                     IER, NREC
      LOGICAL                     FIRST
C     .. Local Arrays ..
      CHARACTER*80                REC(2)
C     .. External Functions ..
      DOUBLE PRECISION            X02AHF, X02AJF, X02AMF
      INTEGER                     P01ABF
      EXTERNAL                    X02AHF, X02AJF, X02AMF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                   ABS, DIMAG, DCMPLX, COS, EXP, LOG,
     *                            MIN, DBLE, SIGN, SIN, SQRT
C     .. Save statement ..
      SAVE                        SAFE, LNSAFE, SAFSIN, RTSAFS, FIRST
C     .. Data statements ..
      DATA                        FIRST/.TRUE./
C     .. Executable Statements ..
      IF (FIRST) THEN
         FIRST = .FALSE.
         SAFE = ONE/X02AMF()
         LNSAFE = LOG(SAFE)
         RECEPS = ONE/X02AJF()
         SAFSIN = MIN(X02AHF(ONE),RECEPS)
         IF (SAFSIN.LT.RECEPS**0.75D0) THEN
C         Assume that SAFSIN is approximately sqrt(RECEPS), in which
C         case IFAIL=4 cannot occur.
            RTSAFS = SAFSIN
         ELSE
C         Set RTSAFS to the argument above which SINE and COSINE will
C         return results of less than half precision, assuming that
C         SAFSIN is approximately equal to RECEPS.
            RTSAFS = SQRT(SAFSIN)
         END IF
      END IF
      NREC = 0
      IER = 0
      X = DBLE(Z)
      Y = DIMAG(Z)
      IF (ABS(Y).GT.SAFSIN) THEN
         IER = 5
         NREC = 2
         WRITE (REC,FMT=99995) Z
         S01EAF = ZERO
      ELSE
         COSY = COS(Y)
         SINY = SIN(Y)
         IF (X.GT.LNSAFE) THEN
            IF (COSY.EQ.ZERO) THEN
               RESR = ZERO
            ELSE
               XPLNCY = X + LOG(ABS(COSY))
               IF (XPLNCY.GT.LNSAFE) THEN
                  IER = 1
                  RESR = SIGN(SAFE,COSY)
               ELSE
                  RESR = SIGN(EXP(XPLNCY),COSY)
               END IF
            END IF
            IF (SINY.EQ.ZERO) THEN
               RESI = ZERO
            ELSE
               XPLNSY = X + LOG(ABS(SINY))
               IF (XPLNSY.GT.LNSAFE) THEN
                  IER = IER + 2
                  RESI = SIGN(SAFE,SINY)
               ELSE
                  RESI = SIGN(EXP(XPLNSY),SINY)
               END IF
            END IF
         ELSE
            EXPX = EXP(X)
            RESR = EXPX*COSY
            RESI = EXPX*SINY
         END IF
         S01EAF = DCMPLX(RESR,RESI)
         IF (IER.EQ.3) THEN
            NREC = 2
            WRITE (REC,FMT=99997) Z
         ELSE IF (ABS(Y).GT.RTSAFS) THEN
            IER = 4
            NREC = 2
            WRITE (REC,FMT=99996) Z
         ELSE IF (IER.EQ.1) THEN
            NREC = 2
            WRITE (REC,FMT=99999) Z
         ELSE IF (IER.EQ.2) THEN
            NREC = 2
            WRITE (REC,FMT=99998) Z
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** Argument Z causes overflow in real part of result:'
     *       ,/4X,'Z = (',1P,D13.5,',',D13.5,')')
99998 FORMAT (1X,'** Argument Z causes overflow in imaginary part of r',
     *       'esult:',/4X,'Z = (',1P,D13.5,',',D13.5,')')
99997 FORMAT (1X,'** Argument Z causes overflow in both real and imagi',
     *       'nary parts of result:',/4X,'Z = (',1P,D13.5,',',D13.5,')')
99996 FORMAT (1X,'** The imaginary part of argument Z is so large that',
     *       ' the result is',/4X,'accurate to less than half precisio',
     *       'n: Z = (',1P,D13.5,',',D13.5,')')
99995 FORMAT (1X,'** The imaginary part of argument Z is so large that',
     *       ' the result has no',/4X,'precision: Z = (',1P,D13.5,',',
     *       D13.5,')')
      END
