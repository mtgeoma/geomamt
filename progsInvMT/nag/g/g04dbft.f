      SUBROUTINE G04DBF(TYPE,NT,TMEAN,RDF,C,LDC,CLEVEL,CIL,CIU,ISIG,
     *                  IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Computes simultaneous confidence intervals for a set of
C     treatment means given the standard errors of the differences.
C
C     It is intended to be used after G04BBF or G04BCF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G04DBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CLEVEL, RDF
      INTEGER           IFAIL, LDC, NT
      CHARACTER         TYPE
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LDC,NT), CIL(NT*(NT-1)/2), CIU(NT*(NT-1)/2),
     *                  TMEAN(NT)
      INTEGER           ISIG(NT*(NT-1)/2)
C     .. Local Scalars ..
      DOUBLE PRECISION  DIFF, P, Q, SD, T
      INTEGER           I, IERROR, IFAULT, IJ, J, NN, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G01FBF, G01FDF, G01FMF
      INTEGER           P01ABF
      EXTERNAL          G01FBF, G01FDF, G01FMF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, SQRT
C     .. Executable Statements ..
      IERROR = 1
      NREC = 1
      IF (NT.LT.2) THEN
         WRITE (P01REC,FMT=99999) NT
      ELSE IF (LDC.LT.NT) THEN
         WRITE (P01REC,FMT=99998) LDC, NT
      ELSE IF (RDF.LT.1.0D0) THEN
         WRITE (P01REC,FMT=99997) RDF
      ELSE IF (CLEVEL.LE.0.0D0 .OR. CLEVEL.GE.1.0D0) THEN
         WRITE (P01REC,FMT=99996) CLEVEL
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.0) THEN
         IF (TYPE.EQ.'T' .OR. TYPE.EQ.'t') THEN
            IFAULT = 1
            Q = SQRT(0.5D0)*G01FMF(CLEVEL,RDF,NT,IFAULT)
            IF (IFAULT.EQ.2) THEN
               IERROR = 3
               WRITE (P01REC,FMT=99993)
               GO TO 60
            END IF
         ELSE IF (TYPE.EQ.'B' .OR. TYPE.EQ.'b') THEN
            NN = NT*(NT-1)/2
            P = (1.0D0-CLEVEL)/DBLE(NN)
            IFAULT = 1
            Q = G01FBF('S',P,RDF,IFAULT)
         ELSE IF (TYPE.EQ.'D' .OR. TYPE.EQ.'d') THEN
            NN = NT*(NT-1)/2
            P = CLEVEL**(1.0D0/DBLE(NN))
            IFAULT = 1
            Q = G01FBF('C',P,RDF,IFAULT)
         ELSE IF (TYPE.EQ.'L' .OR. TYPE.EQ.'l') THEN
            IFAULT = 1
            Q = G01FBF('C',CLEVEL,RDF,IFAULT)
         ELSE IF (TYPE.EQ.'S' .OR. TYPE.EQ.'s') THEN
            IFAULT = 1
            Q = DBLE(NT-1)*G01FDF(CLEVEL,DBLE(NT-1),RDF,IFAULT)
            Q = SQRT(Q)
         ELSE
            IERROR = 1
            WRITE (P01REC,FMT=99995) TYPE
            GO TO 60
         END IF
         IJ = 0
         DO 40 I = 1, NT
            DO 20 J = 1, I - 1
               IJ = IJ + 1
               DIFF = TMEAN(I) - TMEAN(J)
               SD = C(I,J)
               IF (SD.LE.0.0D0) THEN
                  IERROR = 2
                  WRITE (P01REC,FMT=99994) I, J
                  GO TO 60
               END IF
               CIL(IJ) = DIFF - Q*SD
               CIU(IJ) = DIFF + Q*SD
               ISIG(IJ) = 0
               T = ABS(DIFF)/SD
               IF (T.GE.Q) ISIG(IJ) = 1
   20       CONTINUE
   40    CONTINUE
   60    CONTINUE
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, NT .lt. 2: NT = ',I16)
99998 FORMAT (' ** On entry, LDC .lt. NT: LDC = ',I16,' NT = ',I16)
99997 FORMAT (' ** On entry, RDF .lt. 1.0: RDF = ',D13.5)
99996 FORMAT (' ** On entry, CLEVEL .le. 0.0 or .ge. 1.0: CLEVEL = ',
     *       D13.5)
99995 FORMAT (' ** On entry, TYPE is not valid: TYPE = ',A1)
99994 FORMAT (' ** On entry, the ',I16,',',I16,'th element of C .le. 0',
     *       '.0')
99993 FORMAT (' ** Failure in computation of studentized range statist',
     *       'ic')
      END
