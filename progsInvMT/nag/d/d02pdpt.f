      SUBROUTINE D02PDP(IER,SRNAME,NREC,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C$$$$ SUBROUTINE DRKMSG $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C************************************************
C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
C************************************************
C
C  Purpose:      To process error messages and terminate the program
C                in the event of a "catastrophic" failure.
C
C  Input:        IER, SRNAME, NREC
C  Output:       IFAIL
C
C  Common:       Initializes:    none
C                Reads:          /JD02PD/ REC
C                Alters:         none
C
C  Comments:
C  =========
C  The output variable IFAIL is assigned the value of the input variable
C  IER.
C
C  IER = -2  reports a successful call of the subroutine SRNAME and
C            indicates that special action is to be taken elsewhere
C            in the suite.  IFAIL is set to 0 and a return is effected.
C
C  IER = 0   reports a successful call of the subroutine SRNAME.  IFAIL
C            is set and a return is effected.
C
C  0 < IER   a message of NREC records contained in
C            the array REC(*) is written to the standard output channel,
C            OUTCH.  IFAIL is set and a return is effected.
C
C     .. Parameters ..
      LOGICAL           TELL
      PARAMETER         (TELL=.FALSE.)
C     .. Scalar Arguments ..
      INTEGER           IER, IFAIL, NREC
      CHARACTER*(*)     SRNAME
C     .. Arrays in Common ..
      CHARACTER*80      REC(10)
C     .. Local Scalars ..
      INTEGER           I
      LOGICAL           OK
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02PDM
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Common blocks ..
      COMMON            /JD02PD/REC
C     .. Save statement ..
      SAVE              /JD02PD/
C     .. Executable Statements ..
C
C  Check if can continue with integrator.
      OK = (SRNAME.EQ.'D02PDF' .OR. SRNAME.EQ.'D02PCF')
     *     .AND. (IER.EQ.2 .OR. IER.EQ.3 .OR. IER.EQ.4)
C
      IF (IER.GT.0) THEN
         NREC = NREC + 1
         IF (OK) THEN
            WRITE (REC(NREC),FMT='(A)')
     *        ' ** You can continue integrating this problem.'
         ELSE
            WRITE (REC(NREC),FMT='(A)')
     *        ' ** You cannot continue integrating this problem.'
         END IF
      END IF
      DO 20 I = NREC + 1, 10
         REC(I) = ' '
   20 CONTINUE
C
C  TELL D02PDM the status of the routine associated with SRNAME
      CALL D02PDM(TELL,SRNAME,IER)
C
C  Negative IER indicates success so round up to 0
      IER = MAX(0,IER)
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
C
      RETURN
C
      END
