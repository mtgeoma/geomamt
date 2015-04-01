      SUBROUTINE G08RBF(NS,NV,NSUM,Y,IP,X,NX,ICEN,GAMMA,NMAX,TOL,PARVAR,
     *                  NPVAR,IRANK,ZIN,ETA,VAPVEC,PAREST,WORK,LWORK,
     *                  IWA,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08RBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  GAMMA, TOL
      INTEGER           IFAIL, IP, LWORK, NMAX, NPVAR, NS, NSUM, NX
C     .. Array Arguments ..
      DOUBLE PRECISION  ETA(NMAX), PAREST(4*IP+1), PARVAR(NPVAR,IP),
     *                  VAPVEC(NMAX*(NMAX+1)/2), WORK(LWORK), X(NX,IP),
     *                  Y(NSUM), ZIN(NMAX)
      INTEGER           ICEN(NSUM), IRANK(NMAX), IWA(4*NMAX), NV(NS)
C     .. Local Scalars ..
      DOUBLE PRECISION  SS
      INTEGER           I, IERROR, INUM, IOUT, J, LW1, LW2, MMAX, N1,
     *                  N6, NIWA, NPEST, NWA
C     .. Local Arrays ..
      CHARACTER*80      P01REC(4)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G08RBX
C     .. Executable Statements ..
      IERROR = 0
      IF (NS.LT.1) IERROR = 1
      IF (TOL.LE.0.0D0) IERROR = 1
      IF (IP.LT.1) IERROR = 1
      IF (NMAX.LE.IP) IERROR = 1
      IF (NPVAR.LT.IP+1) IERROR = 1
      IF (NX.LT.NSUM) IERROR = 1
      IF (GAMMA.LT.0.0D0) IERROR = 1
      IF (LWORK.LT.NMAX*(IP+1)) IERROR = 1
      IF (IERROR.NE.0) THEN
         WRITE (P01REC,FMT=99999) NS, TOL, NMAX, IP, NPVAR, NX, LWORK,
     *     NSUM, GAMMA
         IOUT = 4
         GO TO 120
      END IF
      N6 = 0
      MMAX = 0
      INUM = 0
      DO 20 I = 1, NS
         N6 = N6 + NV(I)
         IF (NV(I).LT.1) THEN
            INUM = INUM + 1
            IERROR = 1
         END IF
         IF (NV(I).GT.MMAX) MMAX = NV(I)
   20 CONTINUE
      IF (IERROR.NE.0) THEN
         WRITE (P01REC,FMT=99998) INUM
         IOUT = 1
         GO TO 120
      END IF
      IF (NSUM.NE.N6) IERROR = 1
      IF (IERROR.NE.0) THEN
         WRITE (P01REC,FMT=99997) N6, NSUM
         IOUT = 2
         GO TO 120
      END IF
      IF (NMAX.NE.MMAX) IERROR = 1
      IF (IERROR.NE.0) THEN
         WRITE (P01REC,FMT=99996) MMAX, NMAX
         IOUT = 2
         GO TO 120
      END IF
      INUM = 0
      DO 40 I = 1, NSUM
         IF (ICEN(I).NE.0 .AND. ICEN(I).NE.1) THEN
            IERROR = 2
            INUM = INUM + 1
         END IF
   40 CONTINUE
      IF (IERROR.NE.0) THEN
         WRITE (P01REC,FMT=99995) INUM
         IOUT = 1
         GO TO 120
      END IF
C
C     TEST FOR IFAIL = 5
C
      DO 80 J = 1, IP
         SS = X(1,J)
         DO 60 I = 2, NSUM
            IF (X(I,J).NE.SS) GO TO 80
   60    CONTINUE
         IERROR = 5
         WRITE (P01REC,FMT=99994) J, SS
         IOUT = 1
         GO TO 120
   80 CONTINUE
C
      N1 = NMAX*(NMAX+1)/2
      NPEST = 4*IP + 1
      NIWA = 4*NMAX
      NWA = NMAX
      LW1 = 1
      LW2 = LW1 + NMAX*IP
      CALL G08RBX(Y,ICEN,X,NX,GAMMA,NS,NV,NMAX,N1,IP,NSUM,TOL,NPEST,
     *            NPVAR,NWA,NIWA,IRANK,ZIN,ETA,VAPVEC,PAREST,PARVAR,
     *            WORK(LW1),WORK(LW2),IWA,IERROR)
      IF (IERROR.NE.0) GO TO 100
      IFAIL = 0
      GO TO 140
  100 IF (IERROR.EQ.3) WRITE (P01REC,FMT=99993)
     *    '** ON ENTRY, ALL OF THE OBSERVATIONS ARE ADJUDGED TO BE TIED'
      IF (IERROR.EQ.4) WRITE (P01REC,FMT=99992)
     *  '** THE MATRIX X''(B-A)X IS EITHER SINGULAR OR NON POSITIVE DEF'
     *    , 'INITE'
      IOUT = 1
  120 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,IOUT,P01REC)
  140 RETURN
C
99999 FORMAT (' ** ON ENTRY, ONE OR MORE OF THE FOLLOWING PARAMETER VA',
     *  'LUES IS ILLEGAL',/'    NS =',I16,'   TOL =',1P,D16.5,'  NM',
     *  'AX =',I16,/'    IP =',I16,' NPVAR =',I16,'    NX =',I16,/' LW',
     *  'ORK =',I16,'  NSUM =',I16,' GAMMA =',1P,D16.5)
99998 FORMAT (' **',I16,' ELEMENTS OF ARRAY NV ARE LESS THAN OR EQUAL ',
     *  'TO ZERO')
99997 FORMAT (' ** THE SUM OF THE ELEMENTS OF ARRAY NV IS ',I16,/'    ',
     *  'WHICH IS NOT EQUAL TO NSUM, NSUM =',I16)
99996 FORMAT (' ** THE LARGEST SAMPLE SIZE IS ',I16,/'    WHICH IS NOT',
     *  ' EQUAL TO NMAX, NMAX =',I16)
99995 FORMAT (' ** ON ENTRY,',I16,' ELEMENTS OF ARRAY ICEN ARE NOT EQU',
     *  'AL TO 0 OR 1')
99994 FORMAT (' ** ON ENTRY, ALL ELEMENTS IN COLUMN ',I3,' OF X ARE EQ',
     *  'UAL TO ',1P,D12.5)
99993 FORMAT (1X,A,I16,A,A)
99992 FORMAT (1X,A,A)
      END
