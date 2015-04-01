      SUBROUTINE F11GBF(IREVCM,U,V,WORK,LWORK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-----------------------------------------------------------------------
C
C     F11GBF - Solver for the symmetric iterative solver suite
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F11GBF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IREVCM, LWORK
C     .. Array Arguments ..
      DOUBLE PRECISION  U(*), V(*), WORK(LWORK)
C     .. Local Scalars ..
      INTEGER           IFAILL, IFAILM, INFO, INFOL, IREVCX, KILL,
     *                  LWORKL, NREC
      LOGICAL           DONE, FIRST
C     .. Local Arrays ..
      DOUBLE PRECISION  RDATA(20)
      INTEGER           IDATA(20)
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06DBF, F06FBF, F11BAZ, F11GBY, F11GBZ
C     .. Intrinsic Functions ..
      INTRINSIC         SIGN
C     .. Save statement ..
      SAVE              RDATA, IFAILL, IFAILM, INFO, IREVCX, KILL,
     *                  IDATA, LWORKL, DONE, FIRST
C     .. Data statements ..
      DATA              IREVCX, DONE, FIRST/0, .FALSE., .TRUE./
C     .. Executable Statements ..
C
C     Initialize
C
      INFO = 0
      NREC = 0
C
C     First call to F11GBF
C
      IF (FIRST) THEN
         KILL = 1
C
C        Check the input data
C
         IF (IFAIL.NE.0) THEN
            IFAILL = SIGN(1,IFAIL)
         ELSE
            IFAILL = 0
         END IF
         IF (DONE .AND. (IREVCM.EQ.4)) THEN
            INFO = 1
            IREVCM = 4
            IFAILL = IFAILM
            GO TO 20
         ELSE IF (IREVCM.NE.0) THEN
            INFO = -1
         END IF
         IF (INFO.NE.0) THEN
            IREVCM = 4
            GO TO 20
         END IF
C
         DONE = .FALSE.
         CALL F11BAZ(2,IDATA,RDATA,INFO)
         IF (INFO.NE.0) THEN
            INFO = 3
         ELSE IF (LWORK.LT.IDATA(11)) THEN
            INFO = -5
         END IF
         IF (INFO.NE.0) THEN
            IREVCM = 4
            GO TO 20
         END IF
         LWORKL = LWORK
C
C     Subsequent calls to F11GBF
C
      ELSE
C
C        Check the input data
C
         IF (IREVCM.EQ.5) THEN
            KILL = 10
         ELSE IF (IREVCM.EQ.6) THEN
            INFO = 8
            KILL = 100
         ELSE IF (IREVCM.NE.IREVCX) THEN
            INFO = -1
            KILL = 100
         END IF
         IREVCX = KILL*IREVCX
      END IF
C
C     Use Conjugate Gradient
C
      IF (IDATA(2).LE.1) THEN
         CALL F11GBZ(IREVCX,IDATA,RDATA,U,V,WORK,LWORKL,INFO)
C
C     Use Lanczos method (SYMMLQ)
C
      ELSE
         CALL F11GBY(IREVCX,IDATA,RDATA,U,V,WORK,LWORKL,INFO)
      END IF
C
      IREVCM = IREVCX
C
C     Completion
C
   20 CONTINUE
C
C     Prepare for next call
C
      IF (IREVCM.LE.3) THEN
         IFAIL = 0
         IF ((IREVCM.LE.2) .AND. FIRST) FIRST = .FALSE.
         CALL F11BAZ(3,IDATA,RDATA,INFOL)
C
C     Termination
C
      ELSE
         IF (FIRST) THEN
            IF ((INFO.NE.1) .AND. (INFO.NE.3)) THEN
               CALL F11BAZ(2,IDATA,RDATA,INFOL)
               CALL F06DBF(20,0,IDATA,1)
               CALL F06FBF(20,ZERO,RDATA,1)
               CALL F11BAZ(3,IDATA,RDATA,INFOL)
            END IF
         ELSE
            DONE = .TRUE.
            CALL F11BAZ(3,IDATA,RDATA,INFOL)
         END IF
         IF (INFO.NE.0) THEN
            NREC = 2
            IF (INFO.LT.0) THEN
               NREC = 1
               IF (FIRST) THEN
                  WRITE (REC(1),FMT=99999) 'first ', -INFO
               ELSE
                  WRITE (REC(1),FMT=99999) 're-', -INFO
               END IF
            ELSE IF (INFO.EQ.1) THEN
               NREC = 1
               WRITE (REC(1),FMT=99998)
            ELSE IF (INFO.EQ.2) THEN
               IF (KILL.GT.1) THEN
                  WRITE (REC,FMT=99997)
               ELSE
                  WRITE (REC,FMT=99996)
               END IF
            ELSE IF (INFO.EQ.3) THEN
               WRITE (REC,FMT=99995)
            ELSE IF (INFO.EQ.4) THEN
               WRITE (REC,FMT=99994) IDATA(13)
            ELSE IF (INFO.EQ.5) THEN
               NREC = 1
               WRITE (REC(1),FMT=99993) IDATA(13)
            ELSE IF (INFO.EQ.6) THEN
               WRITE (REC,FMT=99992)
            ELSE IF (INFO.EQ.7) THEN
               WRITE (REC,FMT=99991)
            ELSE
               NREC = 1
               WRITE (REC(1),FMT=99990)
            END IF
         END IF
         FIRST = .TRUE.
         IREVCX = 0
         IFAIL = IFAILL
         IFAILM = IFAILL
C
         IFAIL = P01ABF(IFAIL,INFO,SRNAME,NREC,REC)
      END IF
C
C     End of subroutine F11GBF
C
      RETURN
C
99999 FORMAT (' ** On ',A,'entry, parameter number ',I2,' had an illeg',
     *       ' al value.')
99998 FORMAT (' ** F11GBF has already completed its tasks. You need to',
     *       ' set a new problem.')
99997 FORMAT (' ** User-requested termination: the required accuracy c',
     *       'ould not be obtained.',/' ** However, a reasonable accur',
     *       'acy has been achieved.')
99996 FORMAT (' ** The required accuracy could not be obtained.',/' **',
     *       ' However, a reasonable accuracy has been achieved.')
99995 FORMAT (' ** Either F11GAF was not called before calling this ',
     *       'routine',/' ** or it has returned an error.')
99994 FORMAT (' ** User-requested tidy termination.',/' ** The solutio',
     *       'n has not converged after',I7,' iterations.')
99993 FORMAT (' ** The solution has not converged after',I7,' iteratio',
     *       'ns.')
99992 FORMAT (' ** The preconditioner appears not to be positive-defin',
     *       'ite.',/' ** The computation cannot continue.')
99991 FORMAT (' ** The matrix of the coefficients A appears not to be ',
     *       'positive-definite.',/' ** The computation cannot continu',
     *       'e.')
99990 FORMAT (' ** User-requested immediate termination.')
      END
