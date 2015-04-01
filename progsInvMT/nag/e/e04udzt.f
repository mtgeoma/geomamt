      SUBROUTINE E04UDZ(IOPTNS,NOUT,INFORM,OPKEY)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C***********************************************************************
C     E04UDZ  reads the options file from unit  IOPTNS  and loads the
C     options into the relevant elements of the integer and real
C     parameter arrays.
C
C     Systems Optimization Laboratory, Stanford University.
C     This version dated December 18, 1985.
C***********************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           INFORM, IOPTNS, NOUT
C     .. Subroutine Arguments ..
      EXTERNAL          OPKEY
C     .. Local Scalars ..
      INTEGER           IFAIL, NKEY, NREAD
      LOGICAL           PRNT
      CHARACTER*16      KEY
      CHARACTER*72      BUFFER, OLDBUF
C     .. Local Arrays ..
      CHARACTER*16      TOKEN(1)
      CHARACTER*80      REC(9)
C     .. External Subroutines ..
      EXTERNAL          E04UDV, X04BAF, X04BAY, X04BBF
C     .. Executable Statements ..
      PRNT = .TRUE.
C
C     Return if the unit number is out of range.
C
      IF (IOPTNS.LT.0 .OR. IOPTNS.GT.99) THEN
         INFORM = 1
         RETURN
      END IF
C
C     ------------------------------------------------------------------
C     Look for  BEGIN, ENDRUN  or  SKIP.
C     ------------------------------------------------------------------
      NREAD = 0
   20 IFAIL = 1
      CALL X04BBF(IOPTNS,BUFFER,IFAIL)
      IF (IFAIL.NE.0) GO TO 80
C     20 READ (IOPTNS,FMT='(A)',END=80) BUFFER
      NREAD = NREAD + 1
      NKEY = 1
      CALL E04UDV(BUFFER,1,NKEY,TOKEN)
      KEY = TOKEN(1)
      IF (KEY.EQ.'ENDRUN') GO TO 100
      IF (KEY.NE.'BEGIN') THEN
         IF (NREAD.EQ.1 .AND. KEY.NE.'SKIP') THEN
            WRITE (REC,FMT=99999) IOPTNS, BUFFER
            CALL X04BAY(NOUT,9,REC)
         END IF
         GO TO 20
      END IF
C
C     ------------------------------------------------------------------
C     BEGIN found.
C     This is taken to be the first line of an OPTIONS file.
C     Read the second line to see if it is NOLIST.
C     ------------------------------------------------------------------
      OLDBUF = BUFFER
      IFAIL = 1
      CALL X04BBF(IOPTNS,BUFFER,IFAIL)
      IF (IFAIL.NE.0) GO TO 60
C     READ (IOPTNS,FMT='(A)',END=60) BUFFER
C
      CALL OPKEY(NOUT,BUFFER,KEY)
C
      IF (KEY.EQ.'NOLIST') THEN
         PRNT = .FALSE.
      END IF
C
      IF (PRNT) THEN
         WRITE (REC,FMT='(// A / A /)') ' OPTIONS file', ' ------------'
         CALL X04BAY(NOUT,5,REC)
         WRITE (REC,FMT='(6X, A )') OLDBUF, BUFFER
         CALL X04BAY(NOUT,2,REC)
      END IF
C
C     ------------------------------------------------------------------
C     Read the rest of the file.
C     ------------------------------------------------------------------
C     +    while (key .ne. 'end') loop
   40 IF (KEY.NE.'END') THEN
         IFAIL = 1
         CALL X04BBF(IOPTNS,BUFFER,IFAIL)
         IF (IFAIL.NE.0) GO TO 60
C        READ (IOPTNS,FMT='(A)',END=60) BUFFER
         IF (PRNT) THEN
            WRITE (REC,FMT='( 6X, A )') BUFFER
            CALL X04BAF(NOUT,REC(1))
         END IF
C
         CALL OPKEY(NOUT,BUFFER,KEY)
C
         IF (KEY.EQ.'LIST') PRNT = .TRUE.
         IF (KEY.EQ.'NOLIST') PRNT = .FALSE.
         GO TO 40
      END IF
C     +    end while
C
      INFORM = 0
      RETURN
C
   60 WRITE (REC,FMT=99998) IOPTNS
      CALL X04BAY(NOUT,3,REC)
      INFORM = 2
      RETURN
C
   80 WRITE (REC,FMT=99997) IOPTNS
      CALL X04BAY(NOUT,3,REC)
      INFORM = 3
      RETURN
C
  100 WRITE (REC,FMT='(// 6X, A)') BUFFER
      CALL X04BAY(NOUT,3,REC)
      INFORM = 4
      RETURN
C
C
C     End of  E04UDZ. (OPFILE)
C
99999 FORMAT (//' XXX  Error while looking for an OPTIONS file on unit',
     *  I7,/' XXX  The file should start with BEGIN, SKIP or ENDRUN',
     *  /' XXX  but the first record found was the following:',//' ---',
     *  '--',A,//' XXX  Continuing to look for OPTIONS file...')
99998 FORMAT (//' XXX  End-of-file encountered while processing an OPT',
     *  'IONS file on unit',I6)
99997 FORMAT (//' XXX  End-of-file encountered while looking for an OP',
     *  'TIONS file on unit',I6)
      END
