      SUBROUTINE X04CCF(UPLO,DIAG,N,A,TITLE,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     Prints a real triangular matrix stored in packed vector form.
C     Easy to use driver for X04CDF.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='X04CCF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, N
      CHARACTER         DIAG, UPLO
      CHARACTER*(*)     TITLE
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*)
C     .. Local Scalars ..
      INTEGER           IERR
C     .. Local Arrays ..
      CHARACTER         CLABS(1), REC(1), RLABS(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          X04CDF
C     .. Executable Statements ..
      IF (IFAIL.EQ.1) THEN
         IERR = 1
      ELSE
C        Use IERR = -13 so that any error messages from X04CDF do not
C        make reference to X04CDF.
         IERR = -13
      END IF
      CALL X04CDF(UPLO,DIAG,N,A,' ',TITLE,'I',RLABS,'I',CLABS,80,0,IERR)
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,REC)
C
      RETURN
      END
