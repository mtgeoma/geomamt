      SUBROUTINE G04BBF(N,Y,IBLOCK,NT,IT,GMEAN,BMEAN,TMEAN,TABLE,LDT,C,
     *                  LDC,IREP,R,EF,TOL,IRDF,WK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C     MARK 17 REVISED. IER-1662 (JUN 1995).
C
C     Analysis for a block design with one factor.
C     Special cases are:
C                       completely randomised (one-way) design,
C                       randomised complete block,
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G04BBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  GMEAN, TOL
      INTEGER           IBLOCK, IFAIL, IRDF, LDC, LDT, N, NT
C     .. Array Arguments ..
      DOUBLE PRECISION  BMEAN(*), C(LDC,NT), EF(NT), R(N), TABLE(LDT,5),
     *                  TMEAN(NT), WK(3*NT), Y(N)
      INTEGER           IREP(NT), IT(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC, ESS, GAMEAN, RIN, RMS, RSS, SIG, TSS, WKSP
      INTEGER           I, ICREP, IERROR, IFAULT, II, INC, IOBS, IRES,
     *                  ITOTAL, ITREAT, IWARN, J, KBLOCK, NBLOCK, NRDF,
     *                  NREC, NTDF
      LOGICAL           BLOCKS, ORTH
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  G01EDF
      INTEGER           P01ABF
      EXTERNAL          G01EDF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06QHF, G04BBX, G04BBY, G04BBZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, SQRT
C     .. Executable Statements ..
      IERROR = 1
      NREC = 1
      IF (N.LT.2) THEN
         WRITE (P01REC,FMT=99999) N
      ELSE IF (NT.LT.1) THEN
         WRITE (P01REC,FMT=99998) NT
      ELSE IF (NT.LT.2 .AND. ABS(IBLOCK).LE.1) THEN
         WRITE (P01REC,FMT=99997)
      ELSE IF (LDT.LT.4) THEN
         WRITE (P01REC,FMT=99996) LDT
      ELSE IF (LDC.LT.NT) THEN
         WRITE (P01REC,FMT=99995) LDC, NT
      ELSE IF (TOL.LT.0.0D0) THEN
         WRITE (P01REC,FMT=99994) TOL
      ELSE IF (IRDF.LT.0) THEN
         WRITE (P01REC,FMT=99993) IRDF
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.0) THEN
         IF (TOL.EQ.0.0D0) THEN
            ACC = 0.00001D0
         ELSE
            ACC = TOL
         END IF
         ITOTAL = 4
         ITREAT = 2
         IRES = 3
C
C        Inialize TABLE to zero
C
         CALL F06QHF('G',4,5,0.0D0,0.0D0,TABLE,LDT)
C
C        Check blocks
C
         IF (ABS(IBLOCK).GT.1) THEN
            BLOCKS = .TRUE.
            NBLOCK = ABS(IBLOCK)
            KBLOCK = N/NBLOCK
            IF (N.NE.KBLOCK*NBLOCK) THEN
               NREC = 2
               IERROR = 2
               WRITE (P01REC,FMT=99992) N, IBLOCK
               GO TO 420
            END IF
         ELSE
            BLOCKS = .FALSE.
            NBLOCK = 1
            KBLOCK = N
         END IF
C
C        Check treatment factor
C
         IF (NT.GT.1) THEN
            DO 20 J = 1, NT
               IREP(J) = 0
   20       CONTINUE
            ORTH = .TRUE.
            II = 1
            IF (IBLOCK.LT.-1) THEN
               INC = NBLOCK
            ELSE
               INC = 1
            END IF
            DO 100 I = 1, NBLOCK
               DO 40 J = 1, KBLOCK
                  IOBS = IT(II)
                  IF (IOBS.LE.0 .OR. IOBS.GT.NT) THEN
                     IERROR = 3
                     WRITE (P01REC,FMT=99991) II
                     GO TO 420
                  ELSE
                     IREP(IOBS) = IREP(IOBS) + 1
                  END IF
                  II = II + INC
   40          CONTINUE
C
C              Check for equal replication within blocks
C
               IF (ORTH .AND. BLOCKS) THEN
                  ICREP = IREP(1)
                  DO 60 J = 2, NT
                     IF (IREP(J).NE.ICREP) THEN
                        ORTH = .FALSE.
                        GO TO 80
                     END IF
   60             CONTINUE
   80             CONTINUE
               END IF
               IF (INC.NE.1) II = I + 1
  100       CONTINUE
            DO 120 J = 1, NT
               IF (IREP(J).EQ.0) THEN
                  IERROR = 3
                  WRITE (P01REC,FMT=99990) J
                  GO TO 420
               END IF
  120       CONTINUE
         END IF
C
C        Sweep out mean
C
         GMEAN = 0.0D0
         DO 140 I = 1, N
            GMEAN = GMEAN + Y(I)
  140    CONTINUE
         GMEAN = GMEAN/DBLE(N)
         TSS = 0.0D0
         DO 160 I = 1, N
            R(I) = Y(I) - GMEAN
            TSS = TSS + R(I)*R(I)
  160    CONTINUE
         IF (TSS.LE.0.0D0) THEN
            IERROR = 4
            WRITE (P01REC,FMT=99989)
            GO TO 420
         END IF
         TABLE(ITOTAL,2) = TSS
         IF (IRDF.EQ.0) THEN
            TABLE(ITOTAL,1) = N - 1
            NRDF = N - 1
         ELSE
            TABLE(ITOTAL,1) = N - IRDF
            NRDF = N - IRDF
         END IF
         IF (BLOCKS) THEN
C
C           Sweep out blocks if required
C
            CALL G04BBY(N,IBLOCK,KBLOCK,R,BMEAN,RSS,ESS)
            TABLE(1,1) = NBLOCK - 1
            TABLE(1,2) = ESS
            TABLE(1,3) = ESS/DBLE(NBLOCK-1)
            NRDF = NRDF - (NBLOCK-1)
            DO 180 I = 1, NBLOCK
               BMEAN(I) = BMEAN(I) + GMEAN
  180       CONTINUE
         ELSE
            RSS = TSS
         END IF
         IF (NT.GT.1) THEN
            ESS = RSS
C
C           Compute (adjusted) treatment effect means
C
            DO 200 J = 1, NT
               TMEAN(J) = 0.0D0
  200       CONTINUE
            DO 220 I = 1, N
               TMEAN(IT(I)) = TMEAN(IT(I)) + R(I)
  220       CONTINUE
            IF (BLOCKS .AND. ( .NOT. ORTH)) THEN
C
C              Use routine for generalised inverse of reduced matrix
C
               CALL G04BBX(N,IBLOCK,KBLOCK,NT,IT,C,LDC,IREP,EF,NTDF,
     *                     TMEAN,ACC,WK,IWARN)
               IF (IWARN.EQ.-1) THEN
C
C                 An s.e. is zero
C
                  IERROR = 5
                  WRITE (P01REC,FMT=99988)
                  GO TO 420
               ELSE IF (IWARN.EQ.-2) THEN
C
C                 Design disconected
C
                  IERROR = 8
                  WRITE (P01REC,FMT=99987)
               ELSE IF (IWARN.EQ.-3) THEN
C
C                 Treatments totally confounded
C
                  IERROR = 6
                  WRITE (P01REC,FMT=99983)
                  GO TO 420
               ELSE IF (IWARN.EQ.-4) THEN
C
C                 Eigenvalue computation failed to converge
C
                  IERROR = 5
                  WRITE (P01REC,FMT=99986)
                  GO TO 420
               END IF
C
C               Re-compute residuals for non-orthogonal designs
C
               CALL DCOPY(N,Y,1,R,1)
               CALL G04BBZ(N,R,IT,NT,TMEAN,RSS)
               GAMEAN = 0.0D0
               DO 230 I = 1, N
                  GAMEAN = GAMEAN + R(I)
 230              CONTINUE
               GAMEAN = GAMEAN/DBLE(N)
               CALL G04BBY(N,IBLOCK,KBLOCK,R,WK,RSS,WKSP)
            ELSE
C
C              Use simple formulae for one-way and orthogonal designs
C
               DO 240 J = 1, NT
                  TMEAN(J) = TMEAN(J)/DBLE(IREP(J))
  240          CONTINUE
               RIN = 1.0D0/DBLE(N)
               DO 300 I = 1, NT
                  DO 260 J = 1, I - 1
                     C(J,I) = -RIN
  260             CONTINUE
                  C(I,I) = 1.0D0/DBLE(IREP(I)) - RIN
                  DO 280 J = I + 1, NT
                     C(J,I) = SQRT(1.0D0/DBLE(IREP(I))
     *                        +1.0D0/DBLE(IREP(J)))
  280             CONTINUE
  300          CONTINUE
               DO 320 J = 1, NT
                  EF(J) = 1.0D0
  320          CONTINUE
               NTDF = NT - 1
C
C           Sweep out treatment effects
C
               CALL G04BBZ(N,R,IT,NT,TMEAN,RSS)
               GAMEAN = GMEAN
            END IF
C
C           Compute treatment means from effects
C
            DO 340 I = 1, NT
               TMEAN(I) = TMEAN(I) + GAMEAN
  340       CONTINUE
C
C           Compute Treatment SS and df
C
            ESS = ESS - RSS
            IF (ESS.LT.0.0D0) ESS = 0.0D0
            NRDF = NRDF - NTDF
C
C           Compute ANOVA table
C
            TABLE(ITREAT,1) = NTDF
            TABLE(ITREAT,2) = ESS
            TABLE(ITREAT,3) = ESS/DBLE(NTDF)
            TABLE(IRES,1) = NRDF
            TABLE(IRES,2) = RSS
            IF (NRDF.GT.0) THEN
               RMS = RSS/NRDF
               TABLE(IRES,3) = RMS
               IF (RMS.GT.0.0D0) THEN
                  TABLE(ITREAT,4) = TABLE(ITREAT,3)/RMS
                  IFAULT = 1
                  TABLE(ITREAT,5) = G01EDF('U',TABLE(ITREAT,4),
     *                              TABLE(ITREAT,1),DBLE(NRDF),IFAULT)
                  IF (BLOCKS) THEN
                     TABLE(1,4) = TABLE(1,3)/RMS
                     IFAULT = 1
                     TABLE(1,5) = G01EDF('U',TABLE(1,4),TABLE(1,1),
     *                            DBLE(NRDF),IFAULT)
                  END IF
C
C                 Scale se's and generalised inverse
C
                  SIG = SQRT(RMS)
                  DO 400 I = 1, NT
                     DO 360 J = 1, I
                        C(J,I) = C(J,I)*RMS
  360                CONTINUE
                     DO 380 J = I + 1, NT
                        C(J,I) = C(J,I)*SIG
  380                CONTINUE
  400             CONTINUE
               ELSE
                  IERROR = 7
                  WRITE (P01REC,FMT=99985)
               END IF
            ELSE
               IERROR = 7
               WRITE (P01REC,FMT=99984)
            END IF
         END IF
      END IF
  420 CONTINUE
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, N.lt.2: N = ',I16)
99998 FORMAT (' ** On entry, NT.lt.1: NT = ',I16)
99997 FORMAT (' ** On entry, no blocks or treaments in model')
99996 FORMAT (' ** On entry, LDT.lt.4: LDT = ',I16)
99995 FORMAT (' ** On entry, LDC.lt.NT: LDC = ',I16,' NT = ',I16)
99994 FORMAT (' ** On entry, TOL.lt.0.0: TOL = ',D13.5)
99993 FORMAT (' ** On entry, IRDF.lt.0: IRDF = ',I16)
99992 FORMAT (' ** On entry, N not a multiple of number of blocks.',
     *       /'              N = ',I16,' IBLOCK = ',I16)
99991 FORMAT (' ** On entry, the ',I16,'th value of IT is invalid')
99990 FORMAT (' ** On entry, the ',I16,'th treatment is not present')
99989 FORMAT (' ** On entry, the values in Y are constant')
99988 FORMAT (' ** A computed standard error is zero')
99987 FORMAT (' ** The design is disconnected')
99986 FORMAT (' ** The computation of the eigenvalues has failed to co',
     *       'nverge')
99985 FORMAT (' ** The residual mean square is zero')
99984 FORMAT (' ** The residual degrees of freedom is 0')
99983 FORMAT (' ** The treatments are totally confounded with blocks')
      END
