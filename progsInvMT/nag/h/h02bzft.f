      SUBROUTINE H02BZF(N,M,BL,BU,CLAMDA,ISTATE,IWORK,LIWORK,RWORK,
     *                  LRWORK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16 RE-ISSUE. NAG COPYRIGHT 1992.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='H02BZF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LIWORK, LRWORK, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N+M), BU(N+M), CLAMDA(N+M), RWORK(LRWORK)
      INTEGER           ISTATE(N+M), IWORK(LIWORK)
C     .. Local Scalars ..
      INTEGER           IERR, II, IPBU, IPCV, NCTOTL, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
C
      IERR = 0
      NREC = 0
      IF (N.LE.0) THEN
         IERR = 1
         WRITE (REC,FMT=99999) N
         NREC = 2
         GO TO 40
      END IF
      IF (M.LT.0) THEN
         IERR = 1
         WRITE (REC,FMT=99998) M
         NREC = 2
         GO TO 40
      END IF
C
      NCTOTL = N + M
      IPBU = NCTOTL
      IPCV = IPBU + NCTOTL
      DO 20 II = 1, NCTOTL
         ISTATE(II) = IWORK(II)
         BL(II) = RWORK(II)
         BU(II) = RWORK(II+IPBU)
         CLAMDA(II) = RWORK(II+IPCV)
   20 CONTINUE
C
   40 IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** On entry, N.le.0:',/'    N = ',I16)
99998 FORMAT (' ** On entry, M.lt.0:',/'    M = ',I16)
      END
