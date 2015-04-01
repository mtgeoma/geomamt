      SUBROUTINE G08EDF(CL,N,X,M,K,RL,RU,TIL,NGAPS,NCOUNT,EX,CHI,DF,
     *                  PROB,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08EDF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CHI, DF, PROB, RL, RU, TIL
      INTEGER           IFAIL, K, M, N, NGAPS
      CHARACTER*1       CL
C     .. Array Arguments ..
      DOUBLE PRECISION  EX(K), X(N)
      INTEGER           NCOUNT(K)
C     .. Local Scalars ..
      DOUBLE PRECISION  CX, PGAP, SUM
      INTEGER           I, IERROR, IF2, J, LEN, NREC, NSMALL
      LOGICAL           FIRST, INTER, LAST, ONLY
C     .. Local Arrays ..
      CHARACTER*80      P01REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  G01ECF
      INTEGER           P01ABF
      EXTERNAL          G01ECF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Save statement ..
      SAVE              LEN
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
      ELSE IF (N.LT.1) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99998) N
      ELSE IF ((CL.EQ.'S' .OR. CL.EQ.'s') .AND. M.GT.N) THEN
         IERROR = 3
         NREC = 2
         WRITE (P01REC,FMT=99997) M, N
      ELSE IF (K.LE.1) THEN
         IERROR = 4
         WRITE (P01REC,FMT=99996) K
      ELSE IF ((CL.EQ.'S' .OR. CL.EQ.'s') .AND. K.GT.N) THEN
         IERROR = 4
         NREC = 2
         WRITE (P01REC,FMT=99995) K, N
      ELSE IF (RL.GE.RU) THEN
         IERROR = 5
         WRITE (P01REC,FMT=99994) RL, RU
      ELSE IF (TIL.LE.0.0D0) THEN
         IERROR = 5
         WRITE (P01REC,FMT=99993) TIL
      ELSE IF ((RU-RL).GE.TIL) THEN
         IERROR = 5
         NREC = 2
         WRITE (P01REC,FMT=99992) RL, RU, TIL
      ELSE
         NREC = 0
C
C        If this is only call (CL = S or s) or first call (CL = F or f)
C        then initialise NCOUNT and SUM. If this is an intermediate
C        call (CL = I or i) or the final call (CL = L or l) then any
C        gap from the last call is used i.e. LEN is not re-set to 0
C
         IF (ONLY .OR. FIRST) THEN
            NGAPS = 0
            LEN = 1
            DO 20 J = 1, K
               NCOUNT(J) = 0
   20       CONTINUE
         END IF
C
         IF (NGAPS.GE.M .AND. M.GT.0) GO TO 60
         DO 40 I = 1, N
            IF (RL.LT.X(I) .AND. X(I).LT.RU) THEN
               NCOUNT(LEN) = NCOUNT(LEN) + 1
               NGAPS = NGAPS + 1
               LEN = 1
               IF (NGAPS.GE.M .AND. M.GT.0) GO TO 60
            ELSE
               IF (LEN.LE.K-1) LEN = LEN + 1
            END IF
   40    CONTINUE
C
   60    CONTINUE
         IF (ONLY .OR. LAST) THEN
            IF (NGAPS.EQ.0) THEN
               IERROR = 6
               WRITE (P01REC,FMT=99991)
               NREC = 1
               SUM = 0.0D0
            ELSE
               CHI = 0.0D0
               NSMALL = 0
               SUM = DBLE(NGAPS)
               PGAP = (RU-RL)/TIL
               IF (NGAPS.LT.M) THEN
                  IERROR = 8
                  NREC = 1
                  WRITE (P01REC,FMT=99990) NGAPS, M
               END IF
               DO 80 I = 1, K
                  IF (I.LT.K) THEN
                     EX(I) = SUM*PGAP*(1.0D0-PGAP)**(I-1)
                  ELSE
                     EX(I) = SUM*(1.0D0-PGAP)**(I-1)
                  END IF
                  IF (EX(I).LT.1.0D0) NSMALL = NSMALL + 1
                  IF (EX(I).NE.0.0D0) THEN
                     CX = DBLE(NCOUNT(I))
                     CHI = CHI + (CX-EX(I))*(CX-EX(I))/EX(I)
                  ELSE IF (EX(I).EQ.0.0D0 .AND. NCOUNT(I).EQ.0) THEN
                     CHI = CHI + 0.0D0
                  ELSE
                     IERROR = 7
                     NREC = 3
                     WRITE (P01REC,FMT=99989) I
                     GO TO 100
                  END IF
   80          CONTINUE
               IF (NSMALL.NE.0) THEN
                  IERROR = 9
                  NREC = 1
                  WRITE (P01REC,FMT=99988) NSMALL
               END IF
               DF = DBLE(K-1)
               PROB = G01ECF('UPPER',CHI,DF,IF2)
            END IF
         END IF
      END IF
  100 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (1X,'** On entry, CL is not valid : CL = ',A1)
99998 FORMAT (1X,'** On entry, N.lt.1 : N = ',I16)
99997 FORMAT (1X,'** On entry, with CL = ''S'' or ''s'' , M.gt.N :',/4X,
     *       'M = ',I16,'  and N = ',I16)
99996 FORMAT (1X,'** On entry, MAXG.lt.1 : MAXG =',I16)
99995 FORMAT (1X,'** On entry, with CL = ''S'' or ''s'', MAXG.gt.N :',
     *       /4X,'MAXG = ',I16,'  and N = ',I16)
99994 FORMAT (1X,'** On entry, RLO.ge.RUP : RLO = ',D13.5,' and RUP = ',
     *       D13.5)
99993 FORMAT (1X,'** On entry, TOTLEN.le.0.0 : TOTLEN = ',D13.5)
99992 FORMAT (1X,'** On entry, (RUP - RLO).ge.TOTLEN.',/4X,'RLO = ',
     *       D13.5,' RUP = ',D13.5,'  and TOTLEN = ',D13.5)
99991 FORMAT (1X,'** No gaps were found. Larger N or a wider interval ',
     *       'required.')
99990 FORMAT (1X,'** Only ',I16,' gaps were found but ',I16,' were req',
     *       'uested.')
99989 FORMAT (1X,'** The expected frequency in class ',I16,' equals ze',
     *       'ro.',/4X,'The value of (RUP-RLO)/TOTLEN may be too close',
     *       ' to 0.0 or 1.0.',/4X,'or MAXG is too large relative to t',
     *       'he number of gaps found.')
99988 FORMAT (1X,'** ',I16,' classes have expected frequenciesless tha',
     *       'n one.')
      END
