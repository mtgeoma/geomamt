      SUBROUTINE X04DAF(MATRIX,DIAG,M,N,A,LDA,TITLE,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     Prints a general complex matrix.
C     Easy to use driver for X04DBF.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='X04DAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, M, N
      CHARACTER         DIAG, MATRIX
      CHARACTER*(*)     TITLE
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*)
C     .. Local Scalars ..
      INTEGER           IERR
C     .. Local Arrays ..
      CHARACTER         CLABS(1), REC(1), RLABS(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          X04DBF
C     .. Executable Statements ..
      IF (IFAIL.EQ.1) THEN
         IERR = 1
      ELSE
C        Use IERR = -13 so that any error messages from X04DBF do not
C        make reference to X04DBF.
         IERR = -13
      END IF
      CALL X04DBF(MATRIX,DIAG,M,N,A,LDA,'Above',' ',TITLE,'I',RLABS,'I',
     *            CLABS,80,0,IERR)
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,REC)
C
      RETURN
      END
