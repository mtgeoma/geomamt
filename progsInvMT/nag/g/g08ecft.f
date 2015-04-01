      SUBROUTINE G08ECF(CL,N,X,MSIZE,NCOUNT,LDC,EX,CHI,DF,P,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08ECF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CHI, DF, EX, P
      INTEGER           IFAIL, LDC, MSIZE, N
      CHARACTER         CL
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N)
      INTEGER           NCOUNT(LDC,LDC,MSIZE)
C     .. Local Scalars ..
      DOUBLE PRECISION  CX, TX1, TX2, TX3, XNC
      INTEGER           I, IERROR, IF2, IS, IWRK1, IWRK2, J, K, L, M,
     *                  NLEFT, NP, NREC, NTRIPS, NUSED
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
      SAVE              IWRK1, IWRK2, NTRIPS, NLEFT
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
      ELSE IF (N.LT.1 .OR. (ONLY .AND. N.LT.3)) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99998) N
      ELSE IF (MSIZE.LE.1) THEN
         IERROR = 3
         WRITE (P01REC,FMT=99997) MSIZE
      ELSE IF (LDC.LT.MSIZE) THEN
         IERROR = 4
         NREC = 2
         WRITE (P01REC,FMT=99996) LDC, MSIZE
      ELSE
         NREC = 0
         XNC = DBLE(MSIZE)*(1.0D0+X02AJF())
C
C        If this is only call (CL = S or s) or first call (CL = F or f)
C        then initialise NCOUNT and NTRIPS. If this is an intermediate
C        call (CL = I or i) or the final call (CL = L or l) then compare
C        the first 2 values with the last 2 values from the last call
C        and continue.
C
         IF (ONLY .OR. FIRST) THEN
            NTRIPS = 0
            NLEFT = 0
            NUSED = 0
            NP = N
            DO 60 L = 1, MSIZE
               DO 40 K = 1, MSIZE
                  DO 20 J = 1, MSIZE
                     NCOUNT(J,K,L) = 0
   20             CONTINUE
   40          CONTINUE
   60       CONTINUE
         ELSE IF (NLEFT.EQ.2) THEN
            IF (X(1).LT.0.0D0) GO TO 180
            IF (X(1).GT.1.0D0) GO TO 200
            J = IWRK2
            K = IWRK1
            L = INT(XNC*X(1)) + 1
            IF (L.GT.MSIZE) L = MSIZE
            NCOUNT(J,K,L) = NCOUNT(J,K,L) + 1
            NTRIPS = NTRIPS + 1
            NUSED = 1
            NP = N - 1
         ELSE IF (NLEFT.EQ.1) THEN
            IF (N.GT.1) THEN
               IF (X(1).LT.0.0D0 .OR. X(2).LT.0.0D0) GO TO 180
               IF (X(1).GT.1.0D0 .OR. X(2).GT.1.0D0) GO TO 200
               J = IWRK1
               K = INT(XNC*X(1)) + 1
               L = INT(XNC*X(2)) + 1
               IF (K.GT.MSIZE) K = MSIZE
               IF (L.GT.MSIZE) L = MSIZE
               NCOUNT(J,K,L) = NCOUNT(J,K,L) + 1
               NTRIPS = NTRIPS + 1
               NUSED = 2
               NP = N - 2
            ELSE
C
C              If only 1 left and N=1 then add the obs to WRK after
C              checking thats its in range. NLEFT now equals 2.
C
               IWRK2 = IWRK1
               TX1 = X(N)
               IF (TX1.LT.0.0D0) GO TO 180
               IF (TX1.GT.1.0D0) GO TO 200
               IWRK1 = INT(XNC*TX1) + 1
               IF (IWRK1.GT.MSIZE) IWRK1 = MSIZE
               NLEFT = 2
               GO TO 100
            END IF
         ELSE
            NP = N
            NUSED = 0
         END IF
C
         M = N - 2
         IS = NUSED + 1
         DO 80 I = IS, M, 3
            TX1 = X(I)
            TX2 = X(I+1)
            TX3 = X(I+2)
            IF (TX1.LT.0.0D0 .OR. TX2.LT.0.0D0 .OR. TX3.LT.0.0D0)
     *          GO TO 180
            IF (TX1.GT.1.0D0 .OR. TX2.GT.1.0D0 .OR. TX3.GT.1.0D0)
     *          GO TO 200
            J = INT(XNC*TX1) + 1
            K = INT(XNC*TX2) + 1
            L = INT(XNC*TX3) + 1
            IF (J.GT.MSIZE) J = MSIZE
            IF (K.GT.MSIZE) K = MSIZE
            IF (L.GT.MSIZE) L = MSIZE
            NCOUNT(J,K,L) = NCOUNT(J,K,L) + 1
   80    CONTINUE
         NTRIPS = NTRIPS + NP/3
         IF (MOD(NP,3).EQ.2) THEN
            TX1 = X(N)
            TX2 = X(N-1)
            IF (TX1.LT.0.0D0 .OR. TX2.LT.0.0D0) GO TO 180
            IF (TX1.GT.1.0D0 .OR. TX2.GT.1.0D0) GO TO 200
            IWRK1 = INT(XNC*TX1) + 1
            IWRK2 = INT(XNC*TX2) + 1
            IF (IWRK1.GT.MSIZE) IWRK1 = MSIZE
            IF (IWRK2.GT.MSIZE) IWRK2 = MSIZE
            NLEFT = 2
         ELSE IF (MOD(NP,3).EQ.1) THEN
            TX1 = X(N)
            IF (TX1.LT.0.0D0) GO TO 180
            IF (TX1.GT.1.0D0) GO TO 200
            IWRK1 = INT(XNC*TX1) + 1
            IF (IWRK1.GT.MSIZE) IWRK1 = MSIZE
            NLEFT = 1
         ELSE
            NLEFT = 0
         END IF
C
  100    CONTINUE
         IF (ONLY .OR. LAST) THEN
            M = MSIZE*MSIZE*MSIZE
            EX = DBLE(NTRIPS)/DBLE(M)
            IF (EX.EQ.0) THEN
               IERROR = 6
               NREC = 1
               WRITE (P01REC,FMT=99995)
            ELSE
               CHI = 0.0D0
               DO 160 L = 1, MSIZE
                  DO 140 K = 1, MSIZE
                     DO 120 J = 1, MSIZE
                        CX = DBLE(NCOUNT(J,K,L))
                        CHI = CHI + (CX-EX)*(CX-EX)/EX
  120                CONTINUE
  140             CONTINUE
  160          CONTINUE
               DF = DBLE(M-1)
               IF2 = -1
               P = G01ECF('UPPER',CHI,DF,IF2)
               IF (EX.LE.5.0D0) THEN
                  IERROR = 7
                  NREC = 3
                  WRITE (P01REC,FMT=99994) MSIZE, NTRIPS, EX
               END IF
            END IF
         ELSE
            EX = 0.0D0
            CHI = 0.0D0
            DF = 0.0D0
            P = 0.0D0
         END IF
      END IF
      GO TO 220
  180 IERROR = 5
      NREC = 2
      WRITE (P01REC,FMT=99993)
      GO TO 220
  200 IERROR = 5
      NREC = 2
      WRITE (P01REC,FMT=99992)
  220 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (1X,'** On entry, CL is not valid : CL = ',A1)
99998 FORMAT (1X,'** On entry, N.lt.1 or CL = ''S'' or ''s'' and N.lt.',
     *       '3 : N = ',I16)
99997 FORMAT (1X,'** On entry, MSIZE.le.1 : MSIZE = ',I16)
99996 FORMAT (1X,'** On entry, LDC.lt.MSIZE :',/4X,'LDC = ',I16,'  and',
     *       ' MSIZE = ',I16)
99995 FORMAT (1X,'** No triplets were found because .lt. 3 observation',
     *       's were provided.')
99994 FORMAT (1X,'** The expected value for each cell is less than or ',
     *       'equal to 5.0. MSIZE is too',/4X,'large relative to the n',
     *       'umber of triplets : MSIZE =',I16,/4X,'and number of trip',
     *       'lets =',I16,' and expected value =',D13.5)
99993 FORMAT (1X,'** On entry, at least one element of X.lt.0.0 .',/4X,
     *       'All observations must be between 0.0 and 1.0 inclusive.')
99992 FORMAT (1X,'** On entry, at least one element of X.gt.1.0 .',/4X,
     *       'All observations must be between 0.0 and 1.0 inclusive.')
      END
