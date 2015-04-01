      SUBROUTINE G04CAF(N,Y,NFAC,LFAC,NBLOCK,INTER,IRDF,MTERM,TABLE,
     *                  ITOTAL,TMEAN,MAXT,E,IMEAN,SEMEAN,BMEAN,R,IWK,
     *                  IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Computes ANOVA for general complete factorial design, optionally
C     in blocks.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G04CAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, INTER, IRDF, ITOTAL, MAXT, MTERM, N,
     *                  NBLOCK, NFAC
C     .. Array Arguments ..
      DOUBLE PRECISION  BMEAN(NBLOCK+1), E(MAXT), R(N), SEMEAN(MTERM),
     *                  TABLE(MTERM,5), TMEAN(MAXT), Y(N)
      INTEGER           IMEAN(MTERM), IWK(N+3*NFAC), LFAC(NFAC)
C     .. Local Scalars ..
      DOUBLE PRECISION  ESS, RDF, RMS, RSE, RSS, SUM, Y1
      INTEGER           I, IE, IERROR, IFAIL2, IFAULT, II, IR, IREP,
     *                  ITC, ITERM, J, K, KBLOCK, L, NRDF, NREC
      LOGICAL           YCONST
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  G01EDF
      INTEGER           G02EAZ, P01ABF
      EXTERNAL          G01EDF, G02EAZ, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DSCAL, F06QHF, G04CAV, G04CAX, G04CAY, G04CAZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
C
C     Carry out initial error checks
C
      NREC = 1
      IERROR = 1
      IF (N.LT.4) THEN
         WRITE (P01REC,FMT=99999) N
      ELSE IF (NFAC.LT.1) THEN
         WRITE (P01REC,FMT=99998) NFAC
      ELSE IF (NBLOCK.LT.0) THEN
         WRITE (P01REC,FMT=99997) NBLOCK
      ELSE IF (INTER.LT.0) THEN
         WRITE (P01REC,FMT=99996) INTER
      ELSE IF (INTER.GT.NFAC) THEN
         WRITE (P01REC,FMT=99995) INTER, NFAC
      ELSE IF (IRDF.LT.0) THEN
         WRITE (P01REC,FMT=99994) IRDF
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.1) GO TO 280
C
C     Compute number of rows in ANOVA table and number of means
C
      ITOTAL = NFAC
      K = 0
      DO 100 I = 2, INTER
         IFAULT = 0
         ITERM = G02EAZ(I,NFAC,IFAULT)
         ITOTAL = ITOTAL + ITERM
         DO 20 J = 1, I
            IWK(J) = 1
   20    CONTINUE
         DO 40 J = I + 1, NFAC
            IWK(J) = 0
   40    CONTINUE
         DO 80 J = 1, ITERM
            II = 1
            DO 60 L = 1, NFAC
               IF (IWK(L).EQ.1) II = II*LFAC(L)
   60       CONTINUE
            K = K + II
            CALL G04CAV(NFAC,I,IWK,IFAULT)
   80    CONTINUE
  100 CONTINUE
      ITOTAL = ITOTAL + 3
C
C     Check values of MTERM
C
      IF (MTERM.LT.ITOTAL) THEN
         WRITE (P01REC,FMT=99992) MTERM, ITOTAL
         NREC = 2
         IERROR = 2
         GO TO 280
      END IF
C
C     Check values of MAXT
C
      DO 120 I = 1, NFAC
         IF (LFAC(I).LE.1) THEN
            WRITE (P01REC,FMT=99986) I
            IERROR = 2
            GO TO 280
         ELSE
            K = K + LFAC(I)
         END IF
  120 CONTINUE
      IF (MAXT.LT.K) THEN
         NREC = 2
         WRITE (P01REC,FMT=99993) MAXT, K
         IERROR = 2
         IMEAN(1) = K
         GO TO 280
      END IF
      CALL F06QHF('G',ITOTAL,5,0.0D0,0.0D0,TABLE,MTERM)
C
C     Initialize E and TMEAN to zero
C
      DO 140 I = 1, K
         E(I) = 0.0D0
         TMEAN(I) = 0.0D0
  140 CONTINUE
      IF (NBLOCK.GT.1) THEN
C
C        KBLOCK = number of plots per block
C
C        Check block size
C
         KBLOCK = N/NBLOCK
         IF (N.NE.KBLOCK*NBLOCK) THEN
            NREC = 2
            WRITE (P01REC,FMT=99991) N, NBLOCK
            IERROR = 2
            GO TO 280
         END IF
      ELSE
         KBLOCK = N
      END IF
C
C     ITC = number of treatment combinations
C     IREP = plot replication
C
C     Check number of plot treament combinations
C
      ITC = 1
      DO 160 I = 1, NFAC
         ITC = ITC*LFAC(I)
  160 CONTINUE
      IREP = KBLOCK/ITC
      IF (KBLOCK.NE.IREP*ITC) THEN
         NREC = 2
         WRITE (P01REC,FMT=99987)
         IERROR = 2
         GO TO 280
      END IF
C
C     Copy Y to R and remove mean
C
      SUM = 0.0D0
      DO 180 I = 1, N
         SUM = SUM + Y(I)
  180 CONTINUE
      SUM = SUM/DBLE(N)
      RSS = 0.0D0
      Y1 = Y(1)
      YCONST = .TRUE.
      DO 200 I = 1, N
         R(I) = Y(I) - SUM
         RSS = RSS + R(I)*R(I)
         YCONST = YCONST .AND. Y(I) .EQ. Y1
  200 CONTINUE
      BMEAN(1) = SUM
      IF (RSS.LE.0.0D0 .OR. YCONST) THEN
         WRITE (P01REC,FMT=99990)
         IERROR = 3
         GO TO 280
      END IF
      TABLE(ITOTAL,2) = RSS
      TABLE(ITOTAL,3) = 0.0D0
      TABLE(ITOTAL,4) = 0.0D0
      TABLE(ITOTAL,5) = 0.0D0
      IF (IRDF.EQ.0) THEN
         NRDF = N - 1
         TABLE(ITOTAL,1) = N - 1
      ELSE
         NRDF = N - IRDF
         TABLE(ITOTAL,1) = N - IRDF
      END IF
C
C     Sweep out blocks
C
      IF (NBLOCK.GT.1) THEN
         CALL G04CAX(N,NBLOCK,KBLOCK,R,Y,BMEAN(2),ESS)
         TABLE(1,1) = NBLOCK - 1
         TABLE(1,2) = ESS
         NRDF = NRDF - NBLOCK + 1
      END IF
      ITERM = 1
      IE = 0
      IR = 0
C
C     Sweep out main effects
C
      CALL G04CAY(NFAC,LFAC,IREP,N,Y,R,E,TMEAN,IE,SEMEAN,IMEAN,IR,IWK,
     *            TABLE,MTERM,ITERM,RSS,NRDF)
C
C     Sweep out required interactions
C
      IF (INTER.GT.1) THEN
         CALL G04CAZ(INTER,NFAC,LFAC,0,NBLOCK,0,ITC,IREP,N,Y,R,E,TMEAN,
     *               IE,SEMEAN,IMEAN,IR,IWK,TABLE,MTERM,ITERM,RSS,NRDF,
     *               IWK(N+1))
      END IF
C
C     Compute residual sum of squares
C
      RSS = 0.0D0
      DO 220 I = 1, N
         RSS = RSS + R(I)*R(I)
  220 CONTINUE
      TABLE(ITOTAL-1,1) = NRDF
      TABLE(ITOTAL-1,2) = RSS
      IF (NRDF.LE.0) THEN
         WRITE (P01REC,FMT=99989)
         IERROR = 4
      ELSE IF (RSS.LE.0.0D0) THEN
         WRITE (P01REC,FMT=99988)
         IERROR = 4
      ELSE
C
C        Compute mean squares
C
         DO 240 I = 1, ITOTAL - 1
            IF (TABLE(I,1).GT.0.0D0) TABLE(I,3) = TABLE(I,2)/TABLE(I,1)
  240    CONTINUE
C
C        Compute F statistics and probabilities
C
         RMS = TABLE(ITOTAL-1,3)
         RDF = TABLE(ITOTAL-1,1)
         TABLE(ITOTAL-1,4) = 0.0D0
         TABLE(ITOTAL-1,5) = 0.0D0
         DO 260 I = 1, ITOTAL - 2
            TABLE(I,4) = TABLE(I,3)/RMS
            IFAIL2 = 1
            TABLE(I,5) = G01EDF('U',TABLE(I,4),TABLE(I,1),RDF,IFAIL2)
  260    CONTINUE
         RSE = SQRT(RMS)
         ITERM = ITERM - 1
         CALL DSCAL(ITERM,RSE,SEMEAN,1)
      END IF
  280 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, N.lt.4: N = ',I16)
99998 FORMAT (' ** On entry, NFAC.lt.1: NFAC = ',I16)
99997 FORMAT (' ** On entry, NBLOCK.lt.0: NBLOCK = ',I16)
99996 FORMAT (' ** On entry, INTER.lt.0: INTER = ',I16)
99995 FORMAT (' ** On entry, INTER.gt.NFAC: INTER = ',I16,' NFAC = ',
     *       I16)
99994 FORMAT (' ** On entry, IRDF.lt.0: IRDF = ',I16)
99993 FORMAT (' ** On entry, MAXT is too small: MAXT = ',I16,/'    min',
     *       'imum value = ',I16)
99992 FORMAT (' ** On entry, MTERM is too small: MTERM = ',I16,/'    m',
     *       'inimum value = ',I16)
99991 FORMAT (' ** On entry, N is not a multiple of NBLOCK.',/' N = ',
     *       I16,' NBLOCK = ',I16)
99990 FORMAT (' ** On entry all values of Y are identical')
99989 FORMAT (' ** There are no degrees of freedom for the residual su',
     *       'm of squares')
99988 FORMAT (' ** The residual sum of squares is zero')
99987 FORMAT (' ** On entry, the number of plots per block is incompat',
     *       'ible',/' with the number of plot factors')
99986 FORMAT (' ** On entry, the ',I16,'th value of LFAC.le.1')
      END
