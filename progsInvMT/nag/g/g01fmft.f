      DOUBLE PRECISION FUNCTION G01FMF(P,V,IR,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     COMPUTES THE QUANTILE FOR A GIVEN PROBABIITY.
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01FMF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 P, V
      INTEGER                          IFAIL, IR
C     .. Local Scalars ..
      DOUBLE PRECISION                 ABSACC, FX, PX, PY, RELACC, X, Y
      INTEGER                          IER, IERROR, IF1, IFAULT, IND,
     *                                 NREC
C     .. Local Arrays ..
      DOUBLE PRECISION                 C(17)
      CHARACTER*80                     P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION                 G01EMF
      INTEGER                          P01ABF
      EXTERNAL                         G01EMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL                         C05AZF, G01FMZ
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS
C     .. Executable Statements ..
C
      NREC = 1
      G01FMF = 0.0D0
      IERROR = 1
      IF (V.LT.1.0D0) THEN
         WRITE (P01REC,FMT=99999) V
      ELSE IF (IR.LT.2) THEN
         WRITE (P01REC,FMT=99998) IR
      ELSE IF (P.LE.0.0D0 .OR. P.GE.1.0D0) THEN
         WRITE (P01REC,FMT=99997) P
      ELSE
         IERROR = 0
         CALL G01FMZ(P,V,IR,X,Y,PX,PY,IERROR)
         IF (IERROR.EQ.2) THEN
            WRITE (P01REC,FMT=99996)
            GO TO 40
         END IF
         ABSACC = 0.000005D0
         IF (ABS(PX-P).LE.ABSACC) THEN
            G01FMF = X
         ELSE IF (ABS(PY-P).LE.ABSACC) THEN
            G01FMF = Y
         ELSE
            IER = 0
            IND = -1
            FX = PX - P
            C(1) = PY - P
            RELACC = 0.0000025D0
            IFAULT = 1
   20       CONTINUE
            CALL C05AZF(X,Y,FX,ABSACC,IER,C,IND,IFAULT)
            IF (IND.GE.2 .AND. IND.LE.4) THEN
               IF1 = 1
               FX = G01EMF(X,V,IR,IF1) - P
               IF (ABS(FX).GT.ABSACC) GO TO 20
            END IF
            G01FMF = X
            IF (IF1.NE.0 .OR. IFAULT.EQ.5) THEN
               NREC = 2
               IERROR = 3
               WRITE (P01REC,FMT=99995)
            END IF
         END IF
      END IF
   40 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, V.lt.1.0 : V = ',D13.5)
99998 FORMAT (' ** On entry, IR.lt.2 : IR = ',I16)
99997 FORMAT (' ** On entry, P.le.0.0  or  P.ge.1.0 : P = ',D13.5)
99996 FORMAT (' ** Unable to find initial estimate')
99995 FORMAT (' ** Warning - There is some doubt as to whether',/'    ',
     *       '          full accuracy has been achieved')
      END
