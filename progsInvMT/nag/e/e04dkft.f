      SUBROUTINE E04DKF(STRING)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     E04DKF  loads the option supplied in  STRING  into the relevant
C     element of  IPRMLS  or  RPRMLS.
C
C     -- Written on 5-June-1986.
C        Janet Welding and Sven Hammarling, NAG Central Office.
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
      EXTERNAL          E04DGS, X02ZAZ, X04BAF, X04BAY
C     .. Common blocks ..
      COMMON            /AX02ZA/WMACH
      COMMON            /CE04DG/NEWOPT
C     .. Save statement ..
      SAVE              /AX02ZA/, /CE04DG/, FIRST, PRNT
C     .. Data statements ..
      DATA              FIRST/.TRUE./
C     .. Executable Statements ..
C
C     If first time in, set  NOUT.
C     NEWOPT  is true first time into  E04DJF  or  E04DKF
C     and just after a call to  E04DGF.
C     PRNT  is set to true whenever  NEWOPT  is true.
C
      IF (FIRST) THEN
         FIRST = .FALSE.
         NEWOPT = .TRUE.
      END IF
      CALL X02ZAZ
      NOUT = WMACH(11)
      BUFFER = STRING
C
C     Call  E04DGS  to decode the option and set the parameter value.
C     If NEWOPT is true, reset PRNT and test specially for NOLIST.
C
      IF (NEWOPT) THEN
         NEWOPT = .FALSE.
         PRNT = .TRUE.
         CALL E04DGS(NOUT,BUFFER,KEY)
C
         IF (KEY.EQ.'NOLIST') THEN
            PRNT = .FALSE.
         ELSE
            WRITE (REC,FMT='( // A / A / )') ' Calls to E04DKF',
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
         CALL E04DGS(NOUT,BUFFER,KEY)
C
         IF (KEY.EQ.'LIST') THEN
            PRNT = .TRUE.
         ELSE IF (KEY.EQ.'NOLIST') THEN
            PRNT = .FALSE.
         END IF
      END IF
C
      RETURN
C
C     End of  E04DKF. ( CGOPTN )
C
      END
