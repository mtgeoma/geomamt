      SUBROUTINE E04UCP(N,NCLIN,NCNLN,NCTOTL,LITOTL,LWTOTL)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-721 (DEC 1989).
C     MARK 16 REVISED. IER-1084 (JUL 1993).
C     MARK 17 REVISED. IER-1604 (JUN 1995).
C
C     ******************************************************************
C     E04UCP   allocates the addresses of the work arrays for E04UCZ and
C     E04NCZ.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version   14-February-1985.
C     This version of  E04UCP  dated 12-July-1986.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LENLS
      PARAMETER         (LENLS=20)
      INTEGER           LENNP
      PARAMETER         (LENNP=35)
C     .. Scalar Arguments ..
      INTEGER           LITOTL, LWTOTL, N, NCLIN, NCNLN, NCTOTL
C     .. Scalars in Common ..
      INTEGER           LDT, LDZY, LENNAM, NCOLT
C     .. Arrays in Common ..
      INTEGER           LOCLS(LENLS), LOCNP(LENNP)
C     .. Local Scalars ..
      INTEGER           LADX, LANORM, LAQP, LBL, LBU, LC1MUL, LCJAC,
     *                  LCJDX, LCMUL, LCS1, LCS2, LDLAM, LDSLK, LDX,
     *                  LENAQP, LENT, LENZY, LFEATL, LGQ, LGQ1, LGRAD,
     *                  LHCTRL, LHFRWD, LIPERM, LKACTV, LKX, LNEEDC,
     *                  LQPADX, LQPDX, LQPGQ, LQPHZ, LQPTOL, LRHO,
     *                  LRLAM, LRPQ, LRPQ0, LSLK, LSLK1, LT, LWRK1,
     *                  LWRK2, LWRK3, LWTINF, LX1, LZY, MINIW, MINW
C     .. Common blocks ..
      COMMON            /AE04NC/LOCLS
      COMMON            /AE04UC/LOCNP
      COMMON            /BE04NB/LENNAM, LDT, NCOLT, LDZY
C     .. Executable Statements ..
C
      MINIW = LITOTL + 1
      MINW = LWTOTL + 1
C
C     Assign array lengths that depend upon the problem dimensions.
C
      IF (NCLIN+NCNLN.EQ.0) THEN
         LENT = 0
         LENZY = 0
      ELSE
         LENT = LDT*NCOLT
         LENZY = LDZY*LDZY
      END IF
C
      IF (NCNLN.EQ.0) THEN
         LENAQP = 0
      ELSE
         LENAQP = (NCLIN+NCNLN)*N
      END IF
C
      LKACTV = MINIW
      LKX = LKACTV + N
      LNEEDC = LKX + N
      LIPERM = LNEEDC + NCNLN
      MINIW = LIPERM + NCTOTL
C
      LHFRWD = MINW
      LHCTRL = LHFRWD + N
      LANORM = LHCTRL + N
      LQPGQ = LANORM + NCLIN + NCNLN
      LGQ = LQPGQ + N
      LRLAM = LGQ + N
      LT = LRLAM + N
      LZY = LT + LENT
      MINW = LZY + LENZY
C
      LOCLS(1) = LKACTV
      LOCLS(2) = LANORM
      LOCLS(8) = LQPGQ
      LOCLS(9) = LGQ
      LOCLS(10) = LRLAM
      LOCLS(11) = LT
      LOCLS(12) = LZY
C
C     Assign the addresses for the workspace arrays used by  E04UCU.
C
      LQPADX = MINW
      LQPDX = LQPADX + NCLIN + NCNLN
      LRPQ = LQPDX + N
      LRPQ0 = LRPQ + N
      LQPHZ = LRPQ0 + N
      LWTINF = LQPHZ + N
      LWRK1 = LWTINF + NCTOTL
      LQPTOL = LWRK1 + NCTOTL
      MINW = LQPTOL + NCTOTL
C
      LOCLS(3) = LQPADX
      LOCLS(4) = LQPDX
      LOCLS(5) = LRPQ
      LOCLS(6) = LRPQ0
      LOCLS(7) = LQPHZ
      LOCLS(13) = LWTINF
      LOCLS(14) = LWRK1
      LOCLS(15) = LQPTOL
C
C     Assign the addresses for arrays used in E04UCZ.
C
      LAQP = MINW
      LADX = LAQP + LENAQP
      LBL = LADX + NCLIN + NCNLN
      LBU = LBL + NCTOTL
      LDX = LBU + NCTOTL
      LGQ1 = LDX + N
      LFEATL = LGQ1 + N
      LX1 = LFEATL + NCTOTL
      LWRK2 = LX1 + N
      MINW = LWRK2 + NCTOTL
C
      LOCNP(1) = LKX
      LOCNP(2) = LIPERM
      LOCNP(3) = LAQP
      LOCNP(4) = LADX
      LOCNP(5) = LBL
      LOCNP(6) = LBU
      LOCNP(7) = LDX
      LOCNP(8) = LGQ1
      LOCNP(10) = LFEATL
      LOCNP(11) = LX1
      LOCNP(12) = LWRK2
C
      LCS1 = MINW
      LCS2 = LCS1 + NCNLN
      LC1MUL = LCS2 + NCNLN
      LCMUL = LC1MUL + NCNLN
      LCJDX = LCMUL + NCNLN
      LDLAM = LCJDX + NCNLN
      LDSLK = LDLAM + NCNLN
      LRHO = LDSLK + NCNLN
      LWRK3 = LRHO + NCNLN
      LSLK1 = LWRK3 + NCNLN
      LSLK = LSLK1 + NCNLN
      MINW = LSLK + NCNLN
C
      LOCNP(13) = LCS1
      LOCNP(14) = LCS2
      LOCNP(15) = LC1MUL
      LOCNP(16) = LCMUL
      LOCNP(17) = LCJDX
      LOCNP(18) = LDLAM
      LOCNP(19) = LDSLK
      LOCNP(20) = LRHO
      LOCNP(21) = LWRK3
      LOCNP(22) = LSLK1
      LOCNP(23) = LSLK
      LOCNP(24) = LNEEDC
C
      LCJAC = MINW
      LGRAD = LCJAC + NCNLN*N
      MINW = LGRAD + N
C
      LOCNP(25) = LHFRWD
      LOCNP(26) = LHCTRL
      LOCNP(27) = LCJAC
      LOCNP(28) = LGRAD
C
      LITOTL = MINIW - 1
      LWTOTL = MINW - 1
C
      RETURN
C
C     End of  E04UCP. (NPLOC)
C
      END
