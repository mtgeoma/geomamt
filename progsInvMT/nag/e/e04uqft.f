      SUBROUTINE E04UQF(IOPTNS,INFORM)
C     MARK 14 RELEASE.  NAG COPYRIGHT 1989.
C
C     ******************************************************************
C     E04UQF  reads the options file from unit  IOPTNS  and loads the
C     options into the relevant elements of  IPRMNP, RPRMNP, IPRMNL
C     and  RPRMNL.
C
C     If  IOPTNS .lt. 0  or  IOPTNS .gt. 99  then no file is read,
C     otherwise the file associated with unit  IOPTNS  is read.
C
C     Output:
C
C         INFORM = 0  if a complete  OPTIONS  file was found
C                     (starting with  BEGIN  and ending with  END);
C                  1  if  IOPTNS .lt. 0  or  IOPTNS .gt. 99;
C                  2  if  BEGIN  was found, but end-of-file
C                     occurred before  END  was found;
C                  3  if end-of-file occurred before  BEGIN  or
C                     ENDRUN  were found;
C                  4  if  ENDRUN  was found before  BEGIN.
C     ******************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           INFORM, IOPTNS
C     .. Scalars in Common ..
      LOGICAL           NEWOPT
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      INTEGER           NOUT
      LOGICAL           FIRST
C     .. External Subroutines ..
      EXTERNAL          E04UDZ, E04UPU, X02ZAZ
C     .. Common blocks ..
      COMMON            /AX02ZA/WMACH
      COMMON            /EE04UC/NEWOPT
C     .. Save statement ..
      SAVE              /EE04UC/, /AX02ZA/, FIRST, NOUT
C     .. Data statements ..
      DATA              FIRST/.TRUE./
C     .. Executable Statements ..
C
C     NEWOPT is true first time into E04UQF or E04URF
C     and just after a call to E04UPF.
C
      IF (FIRST) THEN
         FIRST = .FALSE.
         NEWOPT = .TRUE.
      END IF
      CALL X02ZAZ
      NOUT = WMACH(11)
C
      CALL E04UDZ(IOPTNS,NOUT,INFORM,E04UPU)
C
      RETURN
C
C     End of  E04UQF.  (NLFILE)
C
      END
