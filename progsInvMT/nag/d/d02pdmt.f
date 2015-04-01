      SUBROUTINE D02PDM(ASK,SRNAME,STATE)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C$$$$ SUBROUTINE RKSIT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C************************************************
C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
C************************************************
C
C  Purpose:      To save or enquire about the status of each
C                subprogram in the suite.
C
C  Input:        ASK, SRNAME
C  Input/output: STATE
C
C
C  Comments:
C  =========
C  SRNAME indicates which routine is of interest in the call to D02PDM.
C
C  If ASK=.FALSE., then the value of STATE (which as used is the error
C  flag) for the routine SRNAME is saved internally. This value of STATE
C  is usually positive.  There are two exceptions:
C    1. SRNAME='D02PVF' and STATE=-1 indicates a completely new problem,
C       so all SAVEd states are cleared.
C    2. STATE=-2 is used by some routines in the suite to indicate
C       circumstances which require special action.
C
C  If ASK=.TRUE., then D02PDM first checks to see if there were any
C  catastrophic errors, that is, a SAVEd state has value 1. This should
C  happen only when the user has ignored previous error warnings. If
C  no catastrophic errors are flagged, then STATE returns the saved
C  state value for the routine specified by SRNAME.
C
C     .. Parameters ..
      INTEGER           STATES, MINUS1
      PARAMETER         (STATES=7,MINUS1=-1)
C     .. Scalar Arguments ..
      INTEGER           STATE
      LOGICAL           ASK
      CHARACTER*(*)     SRNAME
C     .. Local Scalars ..
      INTEGER           I, NAME
C     .. Local Arrays ..
      INTEGER           SVSTA(STATES)
C     .. Save statement ..
      SAVE              SVSTA
C     .. Data statements ..
      DATA              SVSTA/STATES*MINUS1/
C     .. Executable Statements ..
C
      IF (SRNAME.EQ.'D02PVF') THEN
         NAME = 1
      ELSE IF (SRNAME.EQ.'D02PCF') THEN
         NAME = 2
      ELSE IF (SRNAME.EQ.'D02PYF') THEN
         NAME = 3
      ELSE IF (SRNAME.EQ.'D02PZF') THEN
         NAME = 4
      ELSE IF (SRNAME.EQ.'D02PDF') THEN
         NAME = 5
      ELSE IF (SRNAME.EQ.'D02PXF') THEN
         NAME = 6
      ELSE IF (SRNAME.EQ.'D02PWF') THEN
         NAME = 7
      ELSE
         NAME = 0
      END IF
C
C  (Re)initialize if D02PVF is telling D02PDM to do so.
      IF ( .NOT. ASK .AND. NAME.EQ.1) THEN
         IF (STATE.EQ.MINUS1) THEN
            DO 20 I = 1, STATES
               SVSTA(I) = MINUS1
   20       CONTINUE
            GO TO 60
         END IF
      END IF
C
C  Check for 1 on exit from a previous call.
      IF (ASK) THEN
         DO 40 I = 1, STATES
            IF (SVSTA(I).EQ.1) THEN
               STATE = 1
               GO TO 60
            END IF
   40    CONTINUE
      END IF
C
      IF (ASK) THEN
         STATE = SVSTA(NAME)
      ELSE
         SVSTA(NAME) = STATE
      END IF
C
   60 CONTINUE
C
      RETURN
      END
