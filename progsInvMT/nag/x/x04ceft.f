      SUBROUTINE X04CEF(M,N,KL,KU,A,LDA,TITLE,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     Prints a real banded matrix stored in packed form.
C     Easy to use driver for X04CFF.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='X04CEF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, KL, KU, LDA, M, N
      CHARACTER*(*)     TITLE
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*)
C     .. Local Scalars ..
      INTEGER           IERR
C     .. Local Arrays ..
      CHARACTER         CLABS(1), REC(1), RLABS(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          X04CFF
C     .. Executable Statements ..
      IF (IFAIL.EQ.1) THEN
         IERR = 1
      ELSE
C        Use IERR = -13 so that any error messages from X04CFF do not
C        make reference to X04CFF.
         IERR = -13
      END IF
      CALL X04CFF(M,N,KL,KU,A,LDA,' ',TITLE,'I',RLABS,'I',CLABS,80,0,
     *            IERR)
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,REC)
C
      RETURN
      END
