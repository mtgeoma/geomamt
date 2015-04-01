      SUBROUTINE H02BFZ(BLO,BUP,CLAMDA,ISTATE,NCTOTL,IWORK,LIWORK,RWORK,
     *                  LRWORK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C     .. Scalar Arguments ..
      INTEGER           LIWORK, LRWORK, NCTOTL
C     .. Array Arguments ..
      DOUBLE PRECISION  BLO(NCTOTL), BUP(NCTOTL), CLAMDA(NCTOTL),
     *                  RWORK(LRWORK)
      INTEGER           ISTATE(NCTOTL), IWORK(LIWORK)
C     .. Local Scalars ..
      INTEGER           II, IPBU, IPCV
C     .. Executable Statements ..
      IPBU = NCTOTL
      IPCV = IPBU + NCTOTL
      DO 20 II = 1, NCTOTL
         IWORK(II) = ISTATE(II)
         RWORK(II) = BLO(II)
         RWORK(II+IPBU) = BUP(II)
         RWORK(II+IPCV) = CLAMDA(II)
   20 CONTINUE
C
      RETURN
      END
