      SUBROUTINE G04BCF(NREP,NROW,NCOL,Y,NT,IT,GMEAN,TMEAN,TABLE,LDT,C,
     *                  LDC,IREP,RPMEAN,RMEAN,CMEAN,R,EF,TOL,IRDF,WK,
     *                  IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Analysis for a row and column design with one factor.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G04BCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  GMEAN, TOL
      INTEGER           IFAIL, IRDF, LDC, LDT, NCOL, NREP, NROW, NT
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LDC,NT), CMEAN(NREP*NCOL), EF(NT),
     *                  R(NREP*NROW*NCOL), RMEAN(NREP*NROW),
     *                  RPMEAN(NREP), TABLE(LDT,5), TMEAN(NT), WK(3*NT),
     *                  Y(NREP*NROW*NCOL)
      INTEGER           IREP(NT), IT(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC, CSS, ESS, GAMEAN, REP, RIN, RMS, RSS, SIG,
     *                  SQR2, TSS, WKSP
      INTEGER           I, ICREP, IERROR, IFAULT, II, IOBS, IRES,
     *                  ITOTAL, ITREAT, IWARN, J, K, N, NRDF, NREC, NTDF
      LOGICAL           ORTH
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  G01EDF
      INTEGER           P01ABF
      EXTERNAL          G01EDF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, F06QHF, G04BBY, G04BBZ, G04BCX
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
      IERROR = 1
      NREC = 1
      IF (NREP.LT.1) THEN
         WRITE (P01REC,FMT=99992) NREP
      ELSE IF (NROW.LT.2) THEN
         WRITE (P01REC,FMT=99999) NROW
      ELSE IF (NCOL.LT.2) THEN
         WRITE (P01REC,FMT=99997) NCOL
      ELSE IF (NT.LT.1) THEN
         WRITE (P01REC,FMT=99998) NT
      ELSE IF (LDT.LT.6) THEN
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
         ITOTAL = 6
         ITREAT = 4
         IRES = 5
         N = NREP*NROW*NCOL
C
C        Inialize TABLE to zero
C
         CALL F06QHF('G',6,5,0.0D0,0.0D0,TABLE,LDT)
C
C        Check treatment factor
C
         IF (NT.GT.1) THEN
            DO 20 J = 1, NT
               IREP(J) = 0
   20       CONTINUE
            ORTH = .TRUE.
            II = 1
            DO 140 K = 1, NREP
               DO 100 I = 1, NROW
                  DO 40 J = 1, NCOL
                     IOBS = IT(II)
                     IF (IOBS.LE.0 .OR. IOBS.GT.NT) THEN
                        IERROR = 2
                        WRITE (P01REC,FMT=99991) II
                        GO TO 680
                     ELSE
                        IREP(IOBS) = IREP(IOBS) + 1
                     END IF
                     II = II + 1
   40             CONTINUE
C
C                 Check for equal replication within rows
C
                  IF (ORTH) THEN
                     ICREP = IREP(1)
                     DO 60 J = 2, NT
                        IF (IREP(J).NE.ICREP) THEN
                           ORTH = .FALSE.
                           GO TO 80
                        END IF
   60                CONTINUE
   80                CONTINUE
                  END IF
  100          CONTINUE
               DO 120 J = 1, NT
                  IF (IREP(J).EQ.0) THEN
                     IERROR = 2
                     WRITE (P01REC,FMT=99990) J
                     GO TO 680
                  END IF
  120          CONTINUE
  140       CONTINUE
C
C              Check for equal replication within columns
C
            IF (ORTH) THEN
               DO 160 J = 1, NT
                  IREP(J) = 0
  160          CONTINUE
               DO 260 K = 1, NREP
                  DO 240 I = 1, NCOL
                     II = I
                     DO 180 J = 1, NROW
                        IOBS = IT(II)
                        IREP(IOBS) = IREP(IOBS) + 1
                        II = II + NCOL
  180                CONTINUE
                     ICREP = IREP(1)
                     DO 200 J = 2, NT
                        IF (IREP(J).NE.ICREP) THEN
                           ORTH = .FALSE.
                           GO TO 220
                        END IF
  200                CONTINUE
  220                CONTINUE
  240             CONTINUE
  260          CONTINUE
            END IF
         END IF
C
C        Sweep out mean
C
         GMEAN = 0.0D0
         DO 280 I = 1, N
            GMEAN = GMEAN + Y(I)
  280    CONTINUE
         GMEAN = GMEAN/DBLE(N)
         TSS = 0.0D0
         DO 300 I = 1, N
            R(I) = Y(I) - GMEAN
            TSS = TSS + R(I)*R(I)
  300    CONTINUE
         IF (TSS.LE.0.0D0) THEN
            IERROR = 3
            WRITE (P01REC,FMT=99989)
            GO TO 680
         END IF
         TABLE(ITOTAL,2) = TSS
         IF (IRDF.EQ.0) THEN
            TABLE(ITOTAL,1) = N - 1
            NRDF = N - 1
         ELSE
            TABLE(ITOTAL,1) = N - IRDF
            NRDF = N - IRDF
         END IF
C
C           Sweep out replicates
C
         IF (NREP.GT.1) THEN
            CALL G04BBY(N,NREP,NCOL*NROW,R,RPMEAN,RSS,ESS)
            TABLE(1,1) = NREP - 1
            TABLE(1,2) = ESS
            TABLE(1,3) = ESS/DBLE(NREP-1)
            DO 320 I = 1, NREP
               RPMEAN(I) = RPMEAN(I) + GMEAN
  320       CONTINUE
         END IF
C
C           Sweep out columns
C
         CSS = 0.0D0
         DO 340 J = 1, NREP
            CALL G04BBY(NROW*NCOL,-NCOL,NROW,R((J-1)*NROW*NCOL+1),
     *                  CMEAN((J-1)*NCOL+1),RSS,ESS)
            CSS = CSS + ESS
  340    CONTINUE
         TABLE(3,1) = NREP*(NCOL-1)
         TABLE(3,2) = CSS
         TABLE(3,3) = CSS/DBLE(NREP*(NCOL-1))
         DO 360 I = 1, NREP*NCOL
            CMEAN(I) = CMEAN(I) + GMEAN
  360    CONTINUE
C
C           Sweep out rows
C
         CALL G04BBY(N,NREP*NROW,NCOL,R,RMEAN,RSS,ESS)
         TABLE(2,1) = NREP*(NROW-1)
         TABLE(2,2) = ESS
         TABLE(2,3) = ESS/DBLE(NREP*(NROW-1))
         DO 380 I = 1, NREP*NROW
            RMEAN(I) = RMEAN(I) + GMEAN
  380    CONTINUE
         NRDF = NRDF - NREP*(NROW+NCOL-1) + 1
         IF (NT.GT.1) THEN
            ESS = RSS
C
C           Compute (adjusted) treatment effect means
C
            DO 400 J = 1, NT
               TMEAN(J) = 0.0D0
  400       CONTINUE
            DO 420 I = 1, N
               TMEAN(IT(I)) = TMEAN(IT(I)) + R(I)
  420       CONTINUE
            IF ( .NOT. ORTH) THEN
C
C              Use routine for generalised inverse of reduced matrix
C
               CALL G04BCX(N,NREP,NROW,NCOL,NT,IT,C,LDC,IREP,EF,NTDF,
     *                     TMEAN,ACC,WK,IWARN)
               IF (IWARN.EQ.-1) THEN
C
C                 An s.e. is zero
C
                  IERROR = 4
                  WRITE (P01REC,FMT=99988)
                  GO TO 680
               ELSE IF (IWARN.EQ.-2) THEN
C
C                 Design disconected
C
                  IERROR = 7
                  WRITE (P01REC,FMT=99987)
               ELSE IF (IWARN.EQ.-3) THEN
C
C                 Treatments totally confounded
C
                  IERROR = 5
                  WRITE (P01REC,FMT=99983)
                  GO TO 680
               ELSE IF (IWARN.EQ.-4) THEN
C
C                 Eigenvalue computation failed to converge
C
                  IERROR = 4
                  WRITE (P01REC,FMT=99986)
                  GO TO 680
               END IF
C
C           Re-compute residuals for non-orthogonal designs
C
               CALL DCOPY(N,Y,1,R,1)
               CALL G04BBZ(N,R,IT,NT,TMEAN,RSS)
               GAMEAN = 0.0D0
               DO 440 I = 1, N
                  GAMEAN = GAMEAN + R(I)
  440          CONTINUE
               GAMEAN = GAMEAN/DBLE(N)
               CALL G04BBY(N,NREP,NCOL*NROW,R,WK,RSS,WKSP)
               DO 460 I = 1, NREP
                  CALL G04BBY(NROW*NCOL,-NCOL,NROW,R((I-1)*NROW*NCOL+1),
     *                        WK,RSS,WKSP)
  460          CONTINUE
               CALL G04BBY(N,NREP*NROW,NCOL,R,WK,RSS,WKSP)
            ELSE
C
C              Use simple formulae for orthogonal designs
C
               REP = IREP(1)
               SQR2 = SQRT(2.0D0/REP)
               DO 480 J = 1, NT
                  TMEAN(J) = TMEAN(J)/REP
  480          CONTINUE
               RIN = 1.0D0/DBLE(N)
               DO 540 I = 1, NT
                  DO 500 J = 1, I - 1
                     C(J,I) = -RIN
  500             CONTINUE
                  C(I,I) = 1.0D0/REP - RIN
                  DO 520 J = I + 1, NT
                     C(J,I) = SQR2
  520             CONTINUE
  540          CONTINUE
               DO 560 J = 1, NT
                  EF(J) = 1.0D0
  560          CONTINUE
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
            DO 580 I = 1, NT
               TMEAN(I) = TMEAN(I) + GAMEAN
  580       CONTINUE
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
                  IF (NREP.GT.1) THEN
                     II = 1
                  ELSE
                     II = 2
                  END IF
                  DO 600 I = II, 4
                     TABLE(I,4) = TABLE(I,3)/RMS
                     IFAULT = 1
                     TABLE(I,5) = G01EDF('U',TABLE(I,4),TABLE(I,1),
     *                            DBLE(NRDF),IFAULT)
  600             CONTINUE
C
C                 Scale se's and generalised inverse
C
                  SIG = SQRT(RMS)
                  DO 660 I = 1, NT
                     DO 620 J = 1, I
                        C(J,I) = C(J,I)*RMS
  620                CONTINUE
                     DO 640 J = I + 1, NT
                        C(J,I) = C(J,I)*SIG
  640                CONTINUE
  660             CONTINUE
               ELSE
                  IERROR = 6
                  WRITE (P01REC,FMT=99985)
               END IF
            ELSE
               IERROR = 6
               WRITE (P01REC,FMT=99984)
            END IF
         END IF
      END IF
  680 CONTINUE
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, NROW .lt. 2: NROW = ',I16)
99998 FORMAT (' ** On entry, NT .lt. 1: NT = ',I16)
99997 FORMAT (' ** On entry, NCOL .lt. 2: NCOL = ',I16)
99996 FORMAT (' ** On entry, LDT .lt. 6: LDT = ',I16)
99995 FORMAT (' ** On entry, LDC .lt. NT: LDC = ',I16,' NT = ',I16)
99994 FORMAT (' ** On entry, TOL .lt. 0.0: TOL = ',D13.5)
99993 FORMAT (' ** On entry, IRDF .lt. 0: IRDF = ',I16)
99992 FORMAT (' ** On entry, NREP .lt. 1: NREP = ',I16)
99991 FORMAT (' ** On entry, the ',I16,'th element of IT is invalid')
99990 FORMAT (' ** On entry, the ',I16,'th treatment is not present')
99989 FORMAT (' ** On entry, the elements in Y are constant')
99988 FORMAT (' ** A computed standard error is zero')
99987 FORMAT (' ** The design is disconnected')
99986 FORMAT (' ** The computation of the eigenvalues has failed to co',
     *       'nverge')
99985 FORMAT (' ** The residual mean square is zero')
99984 FORMAT (' ** The residual degrees of freedom is 0')
99983 FORMAT (' ** The treatments are totally confounded with blocks')
      END
