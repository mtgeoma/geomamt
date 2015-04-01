      SUBROUTINE X04CAF(MATRIX,DIAG,M,N,A,LDA,TITLE,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     Prints a general real matrix.
C     Easy to use driver for X04CBF.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='X04CAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, M, N
      CHARACTER         DIAG, MATRIX
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
      EXTERNAL          X04CBF
C     .. Executable Statements ..
      IF (IFAIL.EQ.1) THEN
         IERR = 1
      ELSE
C        Use IERR = -13 so that any error messages from X04CBF do not
C        make reference to X04CBF.
         IERR = -13
      END IF
      CALL X04CBF(MATRIX,DIAG,M,N,A,LDA,' ',TITLE,'I',RLABS,'I',CLABS,
     *            80,0,IERR)
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,REC)
C
      RETURN
      END
