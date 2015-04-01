      SUBROUTINE G08AJF(N1,N2,TAIL,U,P,WRK,LWRK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     The tail probaility is returned via P and corresponds to
C     the TAIL option chosen. This routine calculates the exact
C     probability routine using G08AJZ for the case of no ties.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08AJF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  P, U
      INTEGER           IFAIL, LWRK, N1, N2
      CHARACTER*1       TAIL
C     .. Array Arguments ..
      DOUBLE PRECISION  WRK(LWRK)
C     .. Local Scalars ..
      INTEGER           IERROR, IV, NM, NREC, NWRK
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G08AJZ
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MIN
C     .. Executable Statements ..
      NREC = 1
      IF (N1.LT.1 .OR. N2.LT.1) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) N1, N2
      ELSE IF (TAIL.NE.'L' .AND. TAIL.NE.'l' .AND. TAIL.NE.'U' .AND.
     *         TAIL.NE.'u' .AND. TAIL.NE.'T' .AND. TAIL.NE.'t') THEN
         IERROR = 2
         WRITE (P01REC,FMT=99998) TAIL
      ELSE IF (U.LT.0.0D0) THEN
         IERROR = 3
         WRITE (P01REC,FMT=99997) U
      ELSE IF (LWRK.LT.(N1*N2)/2+1) THEN
         NWRK = N1*N2/2 + 1
         IERROR = 4
         NREC = 2
         WRITE (P01REC,FMT=99996) LWRK, NWRK
      ELSE
         IERROR = 0
C
C              Exact significance level for small samples.
C
         NM = N1*N2
         IV = INT(U)
         IF (TAIL.EQ.'L' .OR. TAIL.EQ.'l') THEN
            IF (IV.EQ.NM) THEN
               P = 1.0D0
            ELSE IF ((2*IV).LE.NM) THEN
               CALL G08AJZ(N1,N2,IV,P,WRK,LWRK)
            ELSE
               IV = NM - IV - 1
               CALL G08AJZ(N1,N2,IV,P,WRK,LWRK)
               P = 1.0D0 - P
            END IF
         ELSE IF (TAIL.EQ.'U' .OR. TAIL.EQ.'u') THEN
            IF (IV.EQ.0) THEN
               P = 1.0D0
            ELSE IF ((2*IV).GE.NM) THEN
               IV = NM - IV
               CALL G08AJZ(N1,N2,IV,P,WRK,LWRK)
            ELSE
               IV = IV - 1
               CALL G08AJZ(N1,N2,IV,P,WRK,LWRK)
               P = 1.0D0 - P
            END IF
         ELSE IF (TAIL.EQ.'T' .OR. TAIL.EQ.'t') THEN
            IF ((2*IV).LE.NM) THEN
               CALL G08AJZ(N1,N2,IV,P,WRK,LWRK)
            ELSE
               IV = NM - IV
               CALL G08AJZ(N1,N2,IV,P,WRK,LWRK)
            END IF
            P = 2.0D0*P
            P = MIN(P,1.0D0)
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (1X,'** On entry, N1.lt.1 or N2.lt.1, N1 = ',I16,
     *       ' , N2 = ',I16)
99998 FORMAT (1X,'** On entry, TAIL is not valid: TAIL = ',A1)
99997 FORMAT (1X,'** On entry, U.lt.zero : U = ',D13.5)
99996 FORMAT (1X,'** On entry, the workspace provided in WRK is not la',
     *       'rge enough.',/4X,'LWRK was ',I16,' but should be at leas',
     *       't ',I16)
      END
