      SUBROUTINE G08EBF(CL,N,X,MSIZE,LAG,NCOUNT,LDC,EX,CHI,DF,P,WRK,
     *                  IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08EBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CHI, DF, EX, P
      INTEGER           IFAIL, LAG, LDC, MSIZE, N
      CHARACTER         CL
C     .. Array Arguments ..
      DOUBLE PRECISION  WRK(2*LAG), X(N)
      INTEGER           NCOUNT(LDC,MSIZE)
C     .. Local Scalars ..
      DOUBLE PRECISION  CX, TX1, TX2, XNC
      INTEGER           I, IERROR, IF2, J, K, L, LAG2, LIM, M, MM,
     *                  NLEFT, NP, NPAIRS, NREC, NUSED
      LOGICAL           FIRST, INTER, LAST, ONLY
C     .. Local Arrays ..
      CHARACTER*80      P01REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  G01ECF, X02AJF
      INTEGER           P01ABF
      EXTERNAL          G01ECF, X02AJF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MOD, DBLE
C     .. Save statement ..
      SAVE              NPAIRS, NLEFT
C     .. Executable Statements ..
C
      NREC = 1
      ONLY = CL .EQ. 'S' .OR. CL .EQ. 's'
      FIRST = CL .EQ. 'F' .OR. CL .EQ. 'f'
      INTER = CL .EQ. 'I' .OR. CL .EQ. 'i'
      LAST = CL .EQ. 'L' .OR. CL .EQ. 'l'
      IERROR = 0
      IF ( .NOT. ONLY .AND. .NOT. FIRST .AND. .NOT. INTER .AND. .NOT.
     *    LAST) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) CL
      ELSE IF (N.LT.1 .OR. (ONLY .AND. N.LT.2)) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99998) N
      ELSE IF (MSIZE.LE.1) THEN
         IERROR = 3
         WRITE (P01REC,FMT=99997) MSIZE
      ELSE IF (LAG.LE.0 .OR. (ONLY .AND. LAG.GE.N)) THEN
         IERROR = 4
         NREC = 2
         WRITE (P01REC,FMT=99996) LAG, N
      ELSE IF (LDC.LT.MSIZE) THEN
         IERROR = 5
         NREC = 2
         WRITE (P01REC,FMT=99995) LDC, MSIZE
      ELSE
         NREC = 0
         XNC = DBLE(MSIZE)*(1.0D0+X02AJF())
         LAG2 = 2*LAG
C
C        If this is only call (CL = S or s) or fisrt call (CL = F or f)
C        then initialise NCOUNT and SUM. If this is an intermediate
C        call (CL = I or i) or the final call (CL = L or l) then compare
C        the fisrt LAG values with the last LAG values from the last
C        call and continue.
C
C        NP = Keeps track of observations left for pairing after the
C             possible use of some due to left over observations from
C             a previous call.
C        NUSED = no of observations used in any pairing done due to
C                left over observations (before main loop is reached).
C        NLEFT = the no of observations unused at the end of the call.
C
         IF (ONLY .OR. FIRST) THEN
            NPAIRS = 0
            NLEFT = 0
            NUSED = 0
            NP = N
            DO 40 K = 1, MSIZE
               DO 20 J = 1, MSIZE
                  NCOUNT(J,K) = 0
   20          CONTINUE
   40       CONTINUE
         ELSE IF (NLEFT.GT.0) THEN
C
C           If there are not 2*LAG obs. and this is an intermediate call
C           then add the observations to the workspace, shifting old
C           obs. back.
C
            IF (INTER .AND. (N+NLEFT.LT.LAG2)) THEN
               DO 60 I = 1, N
                  WRK(NLEFT+I) = X(I)
   60          CONTINUE
               NLEFT = N + NLEFT
               GO TO 300
            END IF
C
C           If last call and not enough obs. to form more pairs then go
C           and calculate statistics.
C
            IF (LAST .AND. (N+NLEFT.LE.LAG)) GO TO 200
C
C           Find LAG pairs. If this is the last call and there are not
C           lag pairs available we still find what pairs we can.
C
            LIM = LAG
            IF (LAST .AND. (N+NLEFT.LT.LAG2)) LIM = N + NLEFT - LAG
            DO 80 I = 1, LIM
               IF ((NLEFT-I-LAG).GE.0) THEN
                  TX1 = WRK(I)
                  TX2 = WRK(I+LAG)
               ELSE IF ((NLEFT-I).GE.0) THEN
                  TX1 = WRK(I)
                  TX2 = X(I+LAG-NLEFT)
               ELSE
                  TX1 = X(I-NLEFT)
                  TX2 = X(I+LAG-NLEFT)
               END IF
               IF (TX1.LT.0.0D0 .OR. TX2.LT.0.0D0) GO TO 260
               IF (TX1.GT.1.0D0 .OR. TX2.GT.1.0D0) GO TO 280
               J = INT(XNC*TX1) + 1
               K = INT(XNC*TX2) + 1
               IF (J.GT.MSIZE) J = MSIZE
               IF (K.GT.MSIZE) K = MSIZE
               NCOUNT(J,K) = NCOUNT(J,K) + 1
   80       CONTINUE
            NPAIRS = NPAIRS + LIM
            NUSED = LAG2 - NLEFT
            NP = N - NUSED
C
C           If last call and not enough obs. left compute test stats.
C
            IF (LAST .AND. NP.LE.LAG) GO TO 200
         ELSE
            NP = N
            NUSED = 0
         END IF
C
C        Two sections for LAG.eq.1 or LAG.gt.1
C
         IF (LAG.EQ.1) THEN
            M = N - 1
            LIM = NUSED + 1
            DO 100 I = LIM, M, 2
               TX1 = X(I)
               TX2 = X(I+1)
               IF (TX1.LT.0.0D0 .OR. TX2.LT.0.0D0) GO TO 260
               IF (TX1.GT.1.0D0 .OR. TX2.GT.1.0D0) GO TO 280
               J = INT(XNC*TX1) + 1
               K = INT(XNC*TX2) + 1
               IF (J.GT.MSIZE) J = MSIZE
               IF (K.GT.MSIZE) K = MSIZE
               NCOUNT(J,K) = NCOUNT(J,K) + 1
  100       CONTINUE
            NPAIRS = NPAIRS + NP/2
            IF (MOD(NP,2).NE.0) THEN
               NLEFT = 1
               WRK(1) = X(N)
            ELSE
               NLEFT = 0
            END IF
         ELSE
C
C           First find out how many groups of 2*LAG we have in the
C           observations left.
C
            MM = INT(DBLE(NP+0.25D0)/DBLE(LAG2))
            M = LAG2*(MM-1) + NUSED
            DO 140 L = NUSED, M, LAG2
               DO 120 I = 1, LAG
                  TX1 = X(L+I)
                  TX2 = X(L+I+LAG)
                  IF (TX1.LT.0.0D0 .OR. TX2.LT.0.0D0) GO TO 260
                  IF (TX1.GT.1.0D0 .OR. TX2.GT.1.0D0) GO TO 280
                  J = INT(XNC*TX1) + 1
                  K = INT(XNC*TX2) + 1
                  IF (J.GT.MSIZE) J = MSIZE
                  IF (K.GT.MSIZE) K = MSIZE
                  NCOUNT(J,K) = NCOUNT(J,K) + 1
  120          CONTINUE
  140       CONTINUE
            NPAIRS = NPAIRS + MM*LAG
C
C           If single or last call then find any remaining pairs,
C           before computing statistics. There will be fewer than
C           LAG pairs if any.
C
            IF (ONLY .OR. LAST) THEN
               NLEFT = NP - LAG2*MM
               IF (NLEFT.GT.LAG) THEN
                  LIM = NLEFT - LAG
                  DO 160 I = LIM, 1, -1
                     TX1 = X(N-I-LAG+1)
                     TX2 = X(N-I+1)
                     IF (TX1.LT.0.0D0 .OR. TX2.LT.0.0D0) GO TO 260
                     IF (TX1.GT.1.0D0 .OR. TX2.GT.1.0D0) GO TO 280
                     J = INT(XNC*TX1) + 1
                     K = INT(XNC*TX2) + 1
                     IF (J.GT.MSIZE) J = MSIZE
                     IF (K.GT.MSIZE) K = MSIZE
                     NCOUNT(J,K) = NCOUNT(J,K) + 1
  160             CONTINUE
                  NPAIRS = NPAIRS + LIM
               END IF
            ELSE
C
C              If there is a call to follow store remaining obs. in
C              workspace. NLEFT could be zero.
C
               NLEFT = NP - LAG2*MM
               DO 180 I = 1, NLEFT
                  WRK(I) = X(N-NLEFT+I)
  180          CONTINUE
            END IF
         END IF
C
C        Calculate the expected value and the chi-squared test statistic
C
  200    CONTINUE
         IF (ONLY .OR. LAST) THEN
            EX = DBLE(NPAIRS)/DBLE(MSIZE*MSIZE)
            IF (EX.EQ.0) THEN
               IERROR = 7
               NREC = 1
               WRITE (P01REC,FMT=99994) LAG
            ELSE
               CHI = 0.0D0
               DO 240 K = 1, MSIZE
                  DO 220 J = 1, MSIZE
                     CX = DBLE(NCOUNT(J,K))
                     CHI = CHI + (CX-EX)*(CX-EX)/EX
  220             CONTINUE
  240          CONTINUE
               DF = DBLE(MSIZE*MSIZE-1)
               IF2 = -1
               P = G01ECF('UPPER',CHI,DF,IF2)
               IF (EX.LE.5.0D0) THEN
                  IERROR = 8
                  NREC = 3
                  WRITE (P01REC,FMT=99993) MSIZE, NPAIRS, EX
               END IF
            END IF
         END IF
      END IF
      GO TO 300
  260 IERROR = 6
      NREC = 2
      WRITE (P01REC,FMT=99992)
      GO TO 300
  280 IERROR = 6
      NREC = 2
      WRITE (P01REC,FMT=99991)
  300 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (1X,'** On entry, CL is not valid : CL = ',A1)
99998 FORMAT (1X,'** On entry, N.lt.1 or CL = ''S'' or ''s'' and N.lt.',
     *       '2 : N = ',I16)
99997 FORMAT (1X,'** On entry, MSIZE.le.1 : MSIZE = ',I16)
99996 FORMAT (1X,'** On entry, LAG.lt.1, or CL = S or s and LAG.ge.N :',
     *       /4X,'LAG = ',I16,' and N = ',I16)
99995 FORMAT (1X,'** On entry, LDC.lt.MSIZE :',/4X,'LDC = ',I16,'  and',
     *       ' MSIZE = ',I16)
99994 FORMAT (1X,'** No pairs were found. LAG must be .ge. total n. : ',
     *       'LAG = ',I16)
99993 FORMAT (1X,'** The expected value for each cell is less than or ',
     *       'equal to 5.0. MSIZE is too',/4X,'large relative to the n',
     *       'umber of pairs. MSIZE = ',I16,/4X,
     *       'and number of pairs = ',I16,'  and expected value = ',
     *       D13.5)
99992 FORMAT (1X,'** On entry at least one element of X is .lt. 0.0 .',
     *       /4X,'All observations must be between 0.0 and 1.0 inclusi',
     *       've.')
99991 FORMAT (1X,'** On entry at least one element of X is .gt. 1.0 .',
     *       /4X,'All observations must be between 0.0 and 1.0 inclusi',
     *       've.')
      END
