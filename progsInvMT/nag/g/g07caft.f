      SUBROUTINE G07CAF(TAIL,EQUAL,NX,NY,XMEAN,YMEAN,XSTD,YSTD,CLEVEL,T,
     *                  DF,PROB,DL,DU,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     This routines performs the two sample t-test for both equal
C     and unequal variances.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G07CAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CLEVEL, DF, DL, DU, PROB, T, XMEAN, XSTD, YMEAN,
     *                  YSTD
      INTEGER           IFAIL, NX, NY
      CHARACTER         EQUAL, TAIL
C     .. Local Scalars ..
      DOUBLE PRECISION  D, DIV, PSTD, SD, SD2, SIG, SV, SX, SY, ZNX, ZNY
      INTEGER           IERROR, IFAULT, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G01EBF, G01FBF
      INTEGER           P01ABF
      EXTERNAL          G01EBF, G01FBF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
C
      NREC = 1
      IERROR = 1
      IF (TAIL.NE.'T' .AND. TAIL.NE.'t' .AND. TAIL.NE.'U' .AND. TAIL.NE.
     *    'u' .AND. TAIL.NE.'L' .AND. TAIL.NE.'l') THEN
         WRITE (P01REC,FMT=99999) TAIL
      ELSE IF (EQUAL.NE.'E' .AND. EQUAL.NE.'e' .AND. EQUAL.NE.'U' .AND.
     *         EQUAL.NE.'u') THEN
         WRITE (P01REC,FMT=99998) EQUAL
      ELSE IF (NX.LT.2) THEN
         WRITE (P01REC,FMT=99997) NX
      ELSE IF (NY.LT.2) THEN
         WRITE (P01REC,FMT=99996) NY
      ELSE IF (XSTD.LE.0.0D0) THEN
         WRITE (P01REC,FMT=99995) XSTD
      ELSE IF (YSTD.LE.0.0D0) THEN
         WRITE (P01REC,FMT=99994) YSTD
      ELSE IF (CLEVEL.LE.0.0D0 .OR. CLEVEL.GE.1.0D0) THEN
         WRITE (P01REC,FMT=99993) CLEVEL
      ELSE
         NREC = 0
         IERROR = 0
         D = XMEAN - YMEAN
         ZNX = DBLE(NX)
         ZNY = DBLE(NY)
         IF (EQUAL.EQ.'E' .OR. EQUAL.EQ.'e') THEN
            DF = DBLE(NX+NY-2)
            PSTD = ((ZNX-1.0D0)*XSTD*XSTD+(ZNY-1.0D0)*YSTD*YSTD)/DF
            SD = SQRT(PSTD*DBLE(NX+NY)/DBLE(NX*NY))
            T = D/SD
         ELSE
            SX = XSTD*XSTD/ZNX
            SY = YSTD*YSTD/ZNY
            SD2 = SX + SY
            SD = SQRT(SD2)
            T = D/SD
            DIV = SX*SX/(ZNX-1.0D0) + SY*SY/(ZNY-1.0D0)
            DF = SD2*SD2/DIV
         END IF
         IFAULT = 0
         IF (TAIL.EQ.'T' .OR. TAIL.EQ.'t') THEN
            PROB = G01EBF('S',T,DF,IFAULT)
         ELSE
            PROB = G01EBF(TAIL,T,DF,IFAULT)
         END IF
         IFAULT = 1
         SIG = G01FBF('C',CLEVEL,DF,IFAULT)
         SV = SIG*SD
         DL = D - SV
         DU = D + SV
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, TAIL is an invalid character : TAIL = ',A1)
99998 FORMAT (' ** On entry, EQUAL is an invalid character : EQUAL = ',
     *       A1)
99997 FORMAT (' ** On entry, NX.lt.2 : NX = ',I16)
99996 FORMAT (' ** On entry, NY.lt.2 : NY = ',I16)
99995 FORMAT (' ** On entry, XSTD.le.0.0 : XSTD = ',D13.5)
99994 FORMAT (' ** On entry, YSTD.le.0.0 : YSTD = ',D13.5)
99993 FORMAT (' ** On entry, CLEVEL.le.0.0 or CLEVEL.ge.1.0 : CLEVEL = '
     *       ,D13.5)
      END
