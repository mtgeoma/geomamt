      SUBROUTINE E04UNV(M,N,LITOTL,LWTOTL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     ******************************************************************
C     E04UNV allocates additional addresses for E04UNF.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version   11-May-1988.
C     This version of  E04UNV dated 12-Jul-94.
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
      INTEGER           LF, LFJAC, LFJDX, LYF, MINIW, MINW
C     .. Common blocks ..
      COMMON            /AE04UP/LOCNL
C     .. Executable Statements ..
C
      MINIW = LITOTL + 1
      MINW = LWTOTL + 1
C
      LF = MINW
      LFJAC = LF + M
      LFJDX = LFJAC + M*N
      LYF = LFJDX + M
      MINW = LYF + M
C
      LOCNL(1) = LF
      LOCNL(2) = LFJAC
      LOCNL(3) = LFJDX
      LOCNL(4) = LYF
C
      LITOTL = MINIW - 1
      LWTOTL = MINW - 1
C
      RETURN
C
C     End of  E04UNV.  (NLLOC)
C
      END
