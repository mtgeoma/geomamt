      SUBROUTINE E04MFN(SUBR,MSG,V,LENV)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04MFN  prints the array V in debug format.
C
C     Original version dated 17-Jul-1987.
C     This version of  E04MFN  dated  31-Jan-1988.
C     ******************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           LENV
      CHARACTER*6       SUBR
      CHARACTER*(*)     MSG
C     .. Array Arguments ..
      DOUBLE PRECISION  V(*)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Scalars ..
      INTEGER           I, II
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
C     .. Executable Statements ..
C
      IF (LENV.LE.0) THEN
         WRITE (REC,FMT=99999) SUBR, MSG
         CALL X04BAY(IPRINT,2,REC)
      ELSE
         WRITE (REC,FMT=99998) SUBR, MSG
         CALL X04BAY(IPRINT,2,REC)
         DO 20 I = 1, LENV, 5
            WRITE (REC,FMT=99997) (V(II),II=I,MIN(I+4,LENV))
            CALL X04BAF(IPRINT,REC(1))
   20    CONTINUE
      END IF
C
      RETURN
C
C     End of  E04MFN. (CMMSG1)
C
99999 FORMAT (/' //',A6,'//  ',A)
99998 FORMAT (/' //',A6,'//  ',A,' ... ')
99997 FORMAT (1P,5D15.5)
      END
