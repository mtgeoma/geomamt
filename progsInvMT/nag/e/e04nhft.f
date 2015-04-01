      SUBROUTINE E04NHF(STRING)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04NHF  loads the option supplied in  STRING  into the relevant
C     element of  IPRMLC  or  RPRMLC.
C     ******************************************************************
C
C     .. Scalar Arguments ..
      CHARACTER*(*)     STRING
C     .. Scalars in Common ..
      LOGICAL           NEWOPT
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      INTEGER           NOUT
      LOGICAL           FIRST, PRNT
      CHARACTER*16      KEY
      CHARACTER*72      BUFFER
C     .. Local Arrays ..
      CHARACTER*80      REC(5)
C     .. External Subroutines ..
      EXTERNAL          E04NFN, X02ZAZ, X04BAF, X04BAY
C     .. Common blocks ..
      COMMON            /AX02ZA/WMACH
      COMMON            /BE04MF/NEWOPT
C     .. Save statement ..
      SAVE              /BE04MF/, /AX02ZA/, FIRST, NOUT, PRNT
C     .. Data statements ..
      DATA              FIRST/.TRUE./
C     .. Executable Statements ..
C
C     If first time in, set  NOUT.
C     NEWOPT  is true first time into  E04NGF  or  E04NHF
C     and just after a call to the main routine (e.g. E04NFF).
C     PRNT  is set to  true  whenever  NEWOPT  is true.
C
      IF (FIRST) THEN
         FIRST = .FALSE.
         NEWOPT = .TRUE.
         CALL X02ZAZ
         NOUT = WMACH(11)
      END IF
      BUFFER = STRING
C
C     Call  E04NFN   to decode the option and set the parameter value.
C     If  NEWOPT  is true, reset  PRNT  and test specially for NOLIST.
C
      IF (NEWOPT) THEN
         NEWOPT = .FALSE.
         PRNT = .TRUE.
         CALL E04NFN(NOUT,BUFFER,KEY)
C
         IF (KEY.EQ.'NOLIST') THEN
            PRNT = .FALSE.
         ELSE
            WRITE (REC,FMT='(// A / A /)') ' Calls to E04NHF',
     *        ' ---------------'
            CALL X04BAY(NOUT,5,REC)
            WRITE (REC,FMT='( 6X, A )') BUFFER
            CALL X04BAF(NOUT,REC(1))
         END IF
      ELSE
         IF (PRNT) THEN
            WRITE (REC,FMT='( 6X, A )') BUFFER
            CALL X04BAF(NOUT,REC(1))
         END IF
         CALL E04NFN(NOUT,BUFFER,KEY)
C
         IF (KEY.EQ.'LIST') PRNT = .TRUE.
         IF (KEY.EQ.'NOLIST') PRNT = .FALSE.
      END IF
C
      RETURN
C
C     End of  E04NHF.  (QPPRM)
C
      END
