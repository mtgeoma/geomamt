      SUBROUTINE H02BBP(BLO,BUP,CLAMDA,ISTATE,NCTOTL,IOPTCL,IWORK,
     *                  LIWORK,RWORK,LRWORK)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     .. Scalar Arguments ..
      INTEGER           IOPTCL, LIWORK, LRWORK, NCTOTL
C     .. Array Arguments ..
      DOUBLE PRECISION  BLO(NCTOTL), BUP(NCTOTL), CLAMDA(NCTOTL),
     *                  RWORK(LRWORK)
      INTEGER           ISTATE(NCTOTL,IOPTCL), IWORK(LIWORK)
C     .. Local Scalars ..
      INTEGER           II, IPBU, IPCV
C     .. Executable Statements ..
      IPBU = NCTOTL
      IPCV = IPBU + NCTOTL
      DO 20 II = 1, NCTOTL
         IWORK(II) = ISTATE(II,IOPTCL)
         RWORK(II) = BLO(II)
         RWORK(II+IPBU) = BUP(II)
         RWORK(II+IPCV) = CLAMDA(II)
   20 CONTINUE
C
      RETURN
      END
