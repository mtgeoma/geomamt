      SUBROUTINE E04MFV(CSET,N,NCLIN,LITOTL,LWTOTL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1566 (JUN 1995).
C
C     ******************************************************************
C     E04MFV   allocates the addresses of the work arrays for E04MFZ.
C
C     Note that the arrays ( GQ, CQ ) lie in contiguous areas of
C     workspace.
C
C     Original version written  2-January-1987.
C     This version of  E04MFV  dated  18-Nov-1990.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LENLC
      PARAMETER         (LENLC=20)
C     .. Scalar Arguments ..
      INTEGER           LITOTL, LWTOTL, N, NCLIN
      LOGICAL           CSET
C     .. Scalars in Common ..
      INTEGER           LDQ, LDT, LENNAM, NCOLT
C     .. Arrays in Common ..
      INTEGER           LOCLC(LENLC)
C     .. Local Scalars ..
      INTEGER           LAD, LANORM, LCQ, LD, LENCQ, LENQ, LENRT,
     *                  LFEATU, LGQ, LKACTV, LKX, LQ, LR, LRLAM, LT,
     *                  LWRK, LWTINF, MINIW, MINW
C     .. Common blocks ..
      COMMON            /AE04MF/LOCLC
      COMMON            /BE04NB/LENNAM, LDT, NCOLT, LDQ
C     .. Save statement ..
      SAVE              /AE04MF/
C     .. Executable Statements ..
C
C     ------------------------------------------------------------------
C     Refer to the first free space in the work arrays.
C     ------------------------------------------------------------------
      MINIW = LITOTL + 1
      MINW = LWTOTL + 1
C
C     ------------------------------------------------------------------
C     Integer workspace.
C     ------------------------------------------------------------------
      LKACTV = MINIW
      LKX = LKACTV + N
      MINIW = LKX + N
C
C     ------------------------------------------------------------------
C     Real workspace.
C     Assign array lengths that depend upon the problem dimensions.
C     ------------------------------------------------------------------
      LENRT = LDT*NCOLT
      IF (NCLIN.EQ.0) THEN
         LENQ = 0
      ELSE
         LENQ = LDQ*LDQ
      END IF
C
      IF (CSET) THEN
         LENCQ = N
      ELSE
         LENCQ = 0
      END IF
C
C     ------------------------------------------------------------------
C     We start with arrays that can be preloaded by smart users.
C     ------------------------------------------------------------------
      LFEATU = MINW
      MINW = LFEATU + NCLIN + N
C
C     Next comes stuff used by  E04MFZ  and  E04NFZ.
C
      LANORM = MINW
      LAD = LANORM + NCLIN
      LD = LAD + NCLIN
      LGQ = LD + N
      LCQ = LGQ + N
      LRLAM = LCQ + LENCQ
      LR = LRLAM + N
      LT = LR
      LQ = LT + LENRT
      LWTINF = LQ + LENQ
      LWRK = LWTINF + N + NCLIN
      MINW = LWRK + N + NCLIN
C
C     Load the addresses in LOCLC.
C
      LOCLC(1) = LKACTV
      LOCLC(2) = LKX
C
      LOCLC(3) = LFEATU
      LOCLC(4) = LANORM
      LOCLC(5) = LAD
C
      LOCLC(7) = LD
      LOCLC(8) = LGQ
      LOCLC(9) = LCQ
      LOCLC(10) = LRLAM
      LOCLC(11) = LR
      LOCLC(12) = LT
      LOCLC(13) = LQ
      LOCLC(14) = LWTINF
      LOCLC(15) = LWRK
C
      LITOTL = MINIW - 1
      LWTOTL = MINW - 1
C
      RETURN
C
C     End of  E04MFV.  (LPLOC)
C
      END
