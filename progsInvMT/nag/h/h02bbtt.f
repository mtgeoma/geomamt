      SUBROUTINE H02BBT(M,N,TOLFES,BIGBND,TOLIV,MAXNOD,INTFST,MAXDPT,
     *                  MSGLVL,ITMAX)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 16 RE-ISSUE. NAG COPYRIGHT 1992.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, TOLFES, TOLIV
      INTEGER           INTFST, ITMAX, M, MAXDPT, MAXNOD, MSGLVL, N
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Arrays ..
      CHARACTER*80      REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          X04BAY
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
C     .. Executable Statements ..
      WRITE (REC,FMT=99999)
      CALL X04BAY(NOUT,3,REC)
      IF (INTFST.GT.0) THEN
         WRITE (REC,FMT=99998) M, N, MAXDPT
      ELSE
         WRITE (REC,FMT=99997) M, N, MAXDPT
      END IF
      CALL X04BAY(NOUT,3,REC)
      WRITE (REC,FMT=99996) TOLFES, MSGLVL, BIGBND, X02AJF()
      CALL X04BAY(NOUT,3,REC)
      IF (MAXNOD.GT.0) THEN
         WRITE (REC,FMT=99995) TOLIV, ITMAX, MAXNOD
      ELSE
         WRITE (REC,FMT=99994) TOLIV, ITMAX
      END IF
      CALL X04BAY(NOUT,3,REC)
C
      RETURN
C
C
C
99999 FORMAT (/' Parameters',/' ----------')
99998 FORMAT (/' Linear constraints......',I10,8X,'First integer solut',
     *       'ion..',8X,'ON',/' Variables...............',I10,8X,'Max ',
     *       'depth of the tree...',I10)
99997 FORMAT (/' Linear constraints......',I10,8X,'First integer solut',
     *       'ion..',7X,'OFF',/' Variables...............',I10,8X,'Max',
     *       ' depth of the tree...',I10)
99996 FORMAT (/' Feasibility tolerance...',1P,D10.2,8X,'Print level...',
     *       '..........',I10,/' Infinite bound size.....',1P,D10.2,8X,
     *       'EPS (machine precision).',1P,D10.2)
99995 FORMAT (/' Integer feasibility tol.',1P,D10.2,8X,'Iteration limi',
     *       't.........',I10,/' Max number of nodes.....',I10,/)
99994 FORMAT (/' Integer feasibility tol.',1P,D10.2,8X,'Iteration limi',
     *       't.........',I10,/' Max number of nodes.....',6X,'NONE',/)
      END
