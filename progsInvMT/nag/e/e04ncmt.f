      SUBROUTINE E04NCM(LPROB,N,NCLIN,LITOTL,LWTOTL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-717 (DEC 1989).
C     MARK 16 REVISED. IER-1065 (JUL 1993).
C     MARK 17 REVISED. IER-1577 (JUN 1995).
C
C     ******************************************************************
C     E04NCM   allocates the addresses of the work arrays for  E04NCZ.
C
C     Note that the arrays  ( GQ, CQ )  and  ( RES, RES0, HZ )  lie in
C     contiguous areas of workspace.
C     RES, RES0 and HZ are not needed for LP.
C     CQ is defined when the objective has an explicit linear term.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written  29-October-1984.
C     This version of E04NCM dated 14-May-1986.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LENLS
      PARAMETER         (LENLS=20)
C     .. Scalar Arguments ..
      INTEGER           LITOTL, LPROB, LWTOTL, N, NCLIN
C     .. Scalars in Common ..
      INTEGER           LDT, LDZY, LENNAM, NCOLT
C     .. Arrays in Common ..
      INTEGER           LOCLS(LENLS)
C     .. Local Scalars ..
      INTEGER           LANORM, LAP, LCQ, LENCQ, LENRES, LENT, LENZY,
     *                  LFEATL, LGQ, LHZ, LKACTV, LPX, LRES, LRES0,
     *                  LRLAM, LT, LWRK, LWTINF, LZY, MINIW, MINW
C     .. Common blocks ..
      COMMON            /AE04NC/LOCLS
      COMMON            /BE04NB/LENNAM, LDT, NCOLT, LDZY
C     .. Executable Statements ..
      MINIW = LITOTL + 1
      MINW = LWTOTL + 1
C
C     Assign array lengths that depend upon the problem dimensions.
C
      IF (NCLIN.EQ.0) THEN
         LENT = 0
         LENZY = 0
      ELSE
         LENT = LDT*NCOLT
         LENZY = LDZY*LDZY
      END IF
C
      LENCQ = 0
      IF (LPROB.EQ.2*(LPROB/2)) LENCQ = N
      LENRES = 0
      IF (LPROB.GT.2) LENRES = N
C
      LKACTV = MINIW
      MINIW = LKACTV + N
C
      LANORM = MINW
      LAP = LANORM + NCLIN
      LPX = LAP + NCLIN
      LGQ = LPX + N
      LCQ = LGQ + N
      LRES = LCQ + LENCQ
      LRES0 = LRES + LENRES
      LHZ = LRES0 + LENRES
      LRLAM = LHZ + LENRES
      LT = LRLAM + N
      LZY = LT + LENT
      LWTINF = LZY + LENZY
      LWRK = LWTINF + N + NCLIN
      LFEATL = LWRK + N + NCLIN
      MINW = LFEATL + N + NCLIN
C
      LOCLS(1) = LKACTV
      LOCLS(2) = LANORM
      LOCLS(3) = LAP
      LOCLS(4) = LPX
      LOCLS(5) = LRES
      LOCLS(6) = LRES0
      LOCLS(7) = LHZ
      LOCLS(8) = LGQ
      LOCLS(9) = LCQ
      LOCLS(10) = LRLAM
      LOCLS(11) = LT
      LOCLS(12) = LZY
      LOCLS(13) = LWTINF
      LOCLS(14) = LWRK
      LOCLS(15) = LFEATL
C
      LITOTL = MINIW - 1
      LWTOTL = MINW - 1
C
C     End of  E04NCM. (LSLOC)
C
      RETURN
      END
