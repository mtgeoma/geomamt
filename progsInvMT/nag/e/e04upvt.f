      SUBROUTINE E04UPV(M,N,LITOTL,LWTOTL)
C     MARK 14 RELEASE.  NAG COPYRIGHT 1989.
C
C     ******************************************************************
C     E04UPV allocates additional addresses for NLSUBS.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version   11-May-1988.
C     This version of  E04UPV dated 11-May-1988.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LENNL
      PARAMETER         (LENNL=20)
C     .. Scalar Arguments ..
      INTEGER           LITOTL, LWTOTL, M, N
C     .. Arrays in Common ..
      INTEGER           LOCNL(LENNL)
C     .. Local Scalars ..
      INTEGER           LF, LFJAC, LFJDX, LWRK4, MINIW, MINW
C     .. Common blocks ..
      COMMON            /AE04UP/LOCNL
C     .. Executable Statements ..
C
      MINIW = LITOTL + 1
      MINW = LWTOTL + 1
C
      LF = MINW
      LFJAC = LF + M
      LFJDX = LFJAC + N*M
      LWRK4 = LFJDX + M
      MINW = LWRK4 + M
C
      LOCNL(1) = LF
      LOCNL(2) = LFJAC
      LOCNL(3) = LFJDX
      LOCNL(4) = LWRK4
C
      LITOTL = MINIW - 1
      LWTOTL = MINW - 1
C
      RETURN
C
C     End of  E04UPV.  (NLLOC)
C
      END
