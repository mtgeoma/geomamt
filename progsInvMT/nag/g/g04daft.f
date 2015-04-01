      SUBROUTINE G04DAF(NT,TMEAN,IREP,RMS,RDF,NC,CT,LDCT,EST,TABLE,LDT,
     *                  TOL,USETX,TX,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Computes estimate and sums of squares and F statistics
C     for a set of linear contrasts
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G04DAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RDF, RMS, TOL
      INTEGER           IFAIL, LDCT, LDT, NC, NT
      LOGICAL           USETX
C     .. Array Arguments ..
      DOUBLE PRECISION  CT(LDCT,NC), EST(NC), TABLE(LDT,5), TMEAN(NT),
     *                  TX(NT)
      INTEGER           IREP(NT)
C     .. Local Scalars ..
      DOUBLE PRECISION  E, EPS, EX, SS
      INTEGER           I, IERROR, IFAULT, J, K, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G01EDF, X02AJF
      INTEGER           P01ABF
      EXTERNAL          G01EDF, X02AJF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Executable Statements ..
      IERROR = 1
      NREC = 1
      IF (NC.LT.1) THEN
         WRITE (P01REC,FMT=99999) NC
      ELSE IF (NT.LT.2) THEN
         WRITE (P01REC,FMT=99998) NT
      ELSE IF (LDCT.LT.NT) THEN
         WRITE (P01REC,FMT=99997) LDCT, NT
      ELSE IF (LDT.LT.NC) THEN
         WRITE (P01REC,FMT=99996) LDT, NC
      ELSE IF (RMS.LE.0.0D0) THEN
         WRITE (P01REC,FMT=99995) RMS
      ELSE IF (RDF.LT.1.0D0) THEN
         WRITE (P01REC,FMT=99994) RDF
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.0) THEN
C
C        Check for orthogonality
C
         IF (TOL.LE.0.0D0) THEN
            EPS = X02AJF()
         ELSE
            EPS = TOL
         END IF
         DO 40 I = 1, NC
            E = 0.0D0
            DO 20 J = 1, NT
               E = E + CT(J,I)
   20       CONTINUE
            IF (ABS(E).GT.EPS) THEN
               IERROR = 2
               WRITE (P01REC,FMT=99993) I
            END IF
   40    CONTINUE
         DO 100 I = 1, NC
            DO 80 K = I + 1, NC
               SS = 0.0D0
               DO 60 J = 1, NT
                  SS = SS + CT(J,I)*CT(J,K)
   60          CONTINUE
               IF (SS.GT.EPS) THEN
                  IERROR = 2
                  WRITE (P01REC,FMT=99992) I, K
               END IF
   80       CONTINUE
  100    CONTINUE
         DO 140 I = 1, NC
            SS = 0.0D0
            E = 0.0D0
            EX = 0.0D0
            DO 120 J = 1, NT
               E = E + CT(J,I)*TMEAN(J)
               IF (USETX) EX = EX + CT(J,I)*TX(J)
               SS = SS + CT(J,I)*CT(J,I)/DBLE(IREP(J))
  120       CONTINUE
            IF ( .NOT. USETX) EX = E
            EST(I) = EX
            TABLE(I,1) = 1.0D0
            TABLE(I,2) = E*EX/SS
            TABLE(I,3) = TABLE(I,2)
            TABLE(I,4) = TABLE(I,3)/RMS
            IFAULT = 1
            TABLE(I,5) = G01EDF('U',TABLE(I,4),TABLE(I,1),RDF,IFAULT)
  140    CONTINUE
      END IF
      CONTINUE
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, NC .lt. 1: NC = ',I16)
99998 FORMAT (' ** On entry, NT .lt. 2: NT = ',I16)
99997 FORMAT (' ** On entry, LDCT .lt. NT: LDCT = ',I16,' NT = ',I16)
99996 FORMAT (' ** On entry, LDT .lt. NC: LDT = ',I16,' NC = ',I16)
99995 FORMAT (' ** On entry, RMS .le. 0.0: RMS = ',E13.5)
99994 FORMAT (' ** On entry, RDF .lt. 1.0: RDF = ',E13.5)
99993 FORMAT (' ** The ',I16,'th contrast is not orthogonal to the mean'
     *       )
99992 FORMAT (' ** The ',I16,'th and ',I16,'th contrasts are not ortho',
     *       'gonal')
      END
