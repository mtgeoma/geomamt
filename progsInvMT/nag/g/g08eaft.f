      SUBROUTINE G08EAF(CL,N,X,M,MAXR,NRUNS,NCOUNT,EX,C,LDC,CHI,DF,PROB,
     *                  WRK,LWRK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08EAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CHI, DF, PROB
      INTEGER           IFAIL, LDC, LWRK, M, MAXR, N, NRUNS
      CHARACTER*1       CL
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LDC,MAXR), EX(MAXR), WRK(LWRK), X(N)
      INTEGER           NCOUNT(MAXR)
C     .. Local Scalars ..
      DOUBLE PRECISION  XLAST
      INTEGER           I, IERROR, IF2, IRU, J, LN, LPRUN, NREC, NTOT
      LOGICAL           FIRST, INTER, LAST, ONLY
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  G01ECF
      INTEGER           P01ABF
      EXTERNAL          G01ECF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G08EAY, G08EAZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Save statement ..
      SAVE              IRU, NTOT, XLAST, LPRUN
C     .. Executable Statements ..
C
      NREC = 1
      IERROR = 0
      ONLY = CL .EQ. 'S' .OR. CL .EQ. 's'
      FIRST = CL .EQ. 'F' .OR. CL .EQ. 'f'
      INTER = CL .EQ. 'I' .OR. CL .EQ. 'i'
      LAST = CL .EQ. 'L' .OR. CL .EQ. 'l'
      LN = MAXR*(MAXR+5)/2 + 1
C
      IF ( .NOT. ONLY .AND. .NOT. FIRST .AND. .NOT. INTER .AND. .NOT.
     *    LAST) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) CL
      ELSE IF (N.LT.1 .OR. (ONLY .AND. N.LT.3)) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99998) N
      ELSE IF (M.GT.N .AND. ONLY) THEN
         IERROR = 3
         NREC = 2
         WRITE (P01REC,FMT=99997) M, N
      ELSE IF (MAXR.LT.1) THEN
         IERROR = 4
         WRITE (P01REC,FMT=99996) MAXR
      ELSE IF (MAXR.GE.N .AND. ONLY) THEN
         IERROR = 4
         NREC = 2
         WRITE (P01REC,FMT=99995) MAXR, N
      ELSE IF (LDC.LT.MAXR) THEN
         IERROR = 5
         NREC = 2
         WRITE (P01REC,FMT=99994) LDC, MAXR
      ELSE IF (LWRK.LT.LN) THEN
         IERROR = 6
         NREC = 2
         WRITE (P01REC,FMT=99993) LWRK, LN
      ELSE
         NREC = 0
C
C        If this is only call (CL = 'S' or 's') or first call
C        (CL = 'F' or 'f') then initialise NCOUNT, IRU and NTOT.
C        If this is an intermediate call (CL = 'I' or 'i') or the
C        final call (CL = 'L' or 'l') then compare the first
C        value with the last value from the last call and continue.
C
         IF (ONLY .OR. FIRST) THEN
            IRU = 1
            LPRUN = 1
            NTOT = 0
            NRUNS = 0
            DO 20 I = 1, MAXR
               NCOUNT(I) = 0
   20       CONTINUE
         ELSE
            IF (NRUNS.GE.M .AND. M.GT.0) GO TO 60
            IF (X(1).LT.XLAST) THEN
               NTOT = NTOT + LPRUN
               NCOUNT(IRU) = NCOUNT(IRU) + 1
               IRU = 1
               LPRUN = 1
               NRUNS = NRUNS + 1
               IF (NRUNS.GE.M .AND. M.GT.0) GO TO 60
            ELSE IF (X(1).EQ.XLAST) THEN
               IERROR = 7
               NREC = 1
               WRITE (P01REC,FMT=99992)
               GO TO 100
            ELSE
               LPRUN = LPRUN + 1
               IF (IRU.LT.MAXR) IRU = IRU + 1
            END IF
         END IF
C
C        Count the runs up to length MAXR for this call.
C
         DO 40 J = 2, N
            IF (X(J).LT.X(J-1)) THEN
               NTOT = NTOT + LPRUN
               NCOUNT(IRU) = NCOUNT(IRU) + 1
               IRU = 1
               LPRUN = 1
               NRUNS = NRUNS + 1
               IF (NRUNS.GE.M .AND. M.GT.0) GO TO 60
            ELSE IF (X(J).EQ.X(J-1)) THEN
               IERROR = 7
               NREC = 1
               WRITE (P01REC,FMT=99992)
               GO TO 100
            ELSE
               LPRUN = LPRUN + 1
               IF (IRU.LT.MAXR) IRU = IRU + 1
            END IF
   40    CONTINUE
   60    CONTINUE
         IF (FIRST .OR. INTER) THEN
            XLAST = X(N)
         ELSE
C
C           Compute the expected values and the covariances
C
            IF (NTOT.LT.MAXR) THEN
               IERROR = 8
               NREC = 2
               WRITE (P01REC,FMT=99991) MAXR, NTOT
            ELSE
               IF2 = -1
               CALL G08EAZ(NTOT,MAXR,EX,C,LDC,WRK,LWRK,IF2)
C
C              Compute the test statistic
C
               DO 80 I = 1, MAXR
                  WRK(I) = DBLE(NCOUNT(I)) - EX(I)
   80          CONTINUE
               IF2 = -1
               CALL G08EAY('L',MAXR,C,LDC,WRK,CHI,WRK(MAXR+1),IF2)
               IF (IF2.NE.0) THEN
                  IERROR = 9
                  NREC = 1
                  WRITE (P01REC,FMT=99990)
               ELSE
                  DF = DBLE(MAXR)
                  IF2 = -1
                  PROB = G01ECF('UPPER',CHI,DF,IF2)
                  IF (NRUNS.LT.M) THEN
                     IERROR = 10
                     NREC = 1
                     WRITE (P01REC,FMT=99989) NRUNS, M
                  END IF
               END IF
            END IF
         END IF
      END IF
  100 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (1X,'** On entry, CL is not valid : CL = ',A1)
99998 FORMAT (1X,'** On entry, N.lt.1 or with CL = ''S'' or ''s'' N.lt',
     *       '.3 : N =',I16)
99997 FORMAT (1X,'** On entry, with CL = ''S'' or ''s'', M is greater ',
     *       'than N :',/4X,'M = ',I16,' and N = ',I16)
99996 FORMAT (1X,'** On entry, MAXR.lt.1 : MAXR =',I16)
99995 FORMAT (1X,'** On entry, with CL = ''S'' or ''s'', MAXR is great',
     *       'er than or equal to N',/5X,'MAXR = ',I16,' and N = ',I16)
99994 FORMAT (1X,'** On entry, the leading dimension of COV is less th',
     *       'an MAXR :',/5X,'LDCOV = ',I16,' and MAXR = ',I16)
99993 FORMAT (1X,'** On entry, the size of the work array, defined by ',
     *       'LWRK, is not large enough.',/5X,'LWRK = ',I16,' and shou',
     *       'ld be at least = ',I16)
99992 FORMAT (1X,'** There are ties present in the sample ')
99991 FORMAT (1X,'** The total length of the runs found is less than M',
     *       'AXR',/'    MAXR =',I16,' whereas total length of all run',
     *       's =',I16)
99990 FORMAT (1X,'** COV is not positive definite. Try a smaller value',
     *       ' for MAXR.')
99989 FORMAT (1X,'** Only ',I16,' runs found but ',I16,' requested.')
      END
