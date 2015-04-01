      SUBROUTINE G01ARY(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,POSN,N,W)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Place the character representation of the integer  N
C     right-justified in a field  W  spaces wide starting at
C     position  POSN  in the output line  PLOT .
C
C     The variables  IP, INUM, and  IW  are internal versions
C     of  POSN, N, and  W .  We proceed by extracting the
C     digits of  N, starting with the low-order digit,
C     and stacking them in  DSTK. (  ND  counts the digits.)
C     Once we have collected all the digits (and know that
C     W  spaces are sufficient), we skip over any unneeded
C     spaces, put out a minus sign if needed, and then put out
C     the digits, starting with the high-order one.
C
C     This routine calls G01ARZ   and depends on having digits
C     0 through 9 in consecutive elements of  CHARS (1 to 10).
C     It also assumes that the minus sign is at  CHMIN = 24 in  CHARS.
C
C     .. Scalar Arguments ..
      INTEGER           LDP, MAXPTR, N, OUTPTR, PMAX, POSN, PSUB, W
C     .. Array Arguments ..
      CHARACTER         CHARS(*), PLOT(LDP,*)
C     .. Local Scalars ..
      INTEGER           CH0, CHD, CHMIN, I, INUM, IP, IQ, IW, ND
C     .. Local Arrays ..
      INTEGER           DSTK(20)
C     .. External Subroutines ..
      EXTERNAL          G01ARZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Data statements ..
      DATA              CH0, CHMIN/1, 24/
C     .. Executable Statements ..
C
      IW = W
      IF (N.LT.0) IW = IW - 1
      INUM = ABS(N)
C
C     Extract and stack the digits of  INUM,
C     making sure that  N  fits in  W  spaces.
C
      ND = 1
   20 CONTINUE
      IQ = INUM/10
      DSTK(ND) = INUM - IQ*10
      IF (IQ.NE.0) THEN
         INUM = IQ
         ND = ND + 1
         GO TO 20
      END IF
C
C     Unstack the digits from  DSTK  and put them out.
C     Note that when  N  is negative, a minus sign must be inserted in
C     the space before the first digit.  Decreasing  IW  by 1 in the
C     initialization has provided a space for the minus sign.
C
   40 IP = POSN
      IF (IP.EQ.0) IP = OUTPTR
      IP = IP + IW - ND
      IF (N.LT.0) THEN
         CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,IP,CHMIN)
         IP = IP + 1
      END IF
      DO 60 I = ND, 1, -1
         CHD = CH0 + DSTK(I)
         CALL G01ARZ(PLOT,LDP,PSUB,PMAX,OUTPTR,MAXPTR,CHARS,IP,CHD)
         IP = IP + 1
   60 CONTINUE
   80 RETURN
      END
