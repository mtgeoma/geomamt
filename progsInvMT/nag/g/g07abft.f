      SUBROUTINE G07ABF(N,XMEAN,CLEVEL,TL,TU,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G07ABF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CLEVEL, TL, TU, XMEAN
      INTEGER           IFAIL, N
C     .. Local Scalars ..
      DOUBLE PRECISION  HDF, P, PP, TOL, XN, XN2
      INTEGER           IERROR, IFAULT, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  G01FFF, X02AJF
      INTEGER           P01ABF
      EXTERNAL          G01FFF, X02AJF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, DBLE
C     .. Executable Statements ..
C
      IERROR = 0
      NREC = 0
      IF (N.LE.0) THEN
         IERROR = 1
         NREC = 1
         WRITE (P01REC,FMT=99999) N
      ELSE IF (XMEAN.LT.0.0D0) THEN
         IERROR = 1
         NREC = 1
         WRITE (P01REC,FMT=99998) XMEAN
      ELSE IF (CLEVEL.LE.0.0D0 .OR. CLEVEL.GE.1.0D0) THEN
         IERROR = 1
         NREC = 1
         WRITE (P01REC,FMT=99997) CLEVEL
      ELSE
         XN = DBLE(N)
         XN2 = 2.0D0*XN
         HDF = XN*XMEAN
         P = 0.5D0 - CLEVEL/2.0D0
         TOL = MAX(50.0D0*X02AJF(),0.5D-12)
C
C        Use inverse gamma function
C
         IF (HDF.EQ.0.0D0) THEN
            TL = 0.0D0
         ELSE
            IFAULT = 1
            TL = G01FFF(P,HDF,2.0D0,TOL,IFAULT)/XN2
            IF (IFAULT.EQ.5) THEN
               IERROR = 2
               NREC = 2
               WRITE (P01REC,FMT=99996)
               TU = 0.0D0
               GO TO 20
            END IF
         END IF
C
         HDF = HDF + 1.0D0
         PP = CLEVEL + P
         IFAULT = 1
         TU = G01FFF(PP,HDF,2.0D0,TOL,IFAULT)/XN2
         IF (IFAULT.EQ.5) THEN
            IERROR = 2
            NREC = 2
            WRITE (P01REC,FMT=99996)
            TL = 0.0D0
         END IF
C
      END IF
C
   20 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, N.le.0 : N = ',I16)
99998 FORMAT (' ** On entry, XMEAN.lt.0.0 : XMEAN = ',D13.5)
99997 FORMAT (' ** On entry, CLEVEL.le.0.0 or CLEVEL.ge.1.0 : CLEVEL = '
     *       ,D13.5)
99996 FORMAT (' ** When using the relationship with the gamma distribu',
     *       'tion the series to',/'    calculate the gamma probabilit',
     *       'ies has failed to converge')
      END
