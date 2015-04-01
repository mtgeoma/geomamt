      SUBROUTINE G04CAZ(INTER,NFAC,LFAC,MPF,NBLOCK,MREP,ITC,NREPL,N,Y,R,
     *                  E,TMEAN,IE,SEMEAN,IMEAN,IR,IFAC,TABLE,LDT,ITERM,
     *                  RSS,NRDF,IWK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Sweeps out interactions
C
C     Can be used for different strata with MPF and MREP non-zero
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RSS
      INTEGER           IE, INTER, IR, ITC, ITERM, LDT, MPF, MREP, N,
     *                  NBLOCK, NFAC, NRDF, NREPL
C     .. Array Arguments ..
      DOUBLE PRECISION  E(*), R(N), SEMEAN(*), TABLE(LDT,5), TMEAN(*),
     *                  Y(N)
      INTEGER           IFAC(N), IMEAN(*), IWK(3*NFAC), LFAC(NFAC)
C     .. Local Scalars ..
      DOUBLE PRECISION  ESS
      INTEGER           I, ICOMB, IFAIL2, IFAULT, II, ILEVEL, J, JJ, K,
     *                  KBLOCK, KK, KPLOT, L, LL, NDF, NTCPB, NTERM
C     .. External Functions ..
      INTEGER           G02EAZ
      EXTERNAL          G02EAZ
C     .. External Subroutines ..
      EXTERNAL          G04CAU, G04CAV, G04CAW
C     .. Executable Statements ..
      IF (NBLOCK.GT.0) THEN
         KBLOCK = N/NBLOCK
      ELSE
         KBLOCK = N
      END IF
      KPLOT = ITC*NREPL
C
C        NTCPB = number of times the treatment combinations (and
C                 replicates) are repeated per block
C
      IF (MREP.GT.1) THEN
         NTCPB = KBLOCK/MREP/KPLOT
      ELSE IF (MPF.GT.0) THEN
         NTCPB = KBLOCK/KPLOT
      ELSE
         NTCPB = 1
      END IF
C
C        Initialise IWK for (SIMDO)
C
      IWK(2*NFAC+1) = LFAC(NFAC)
      DO 20 I = 2, NFAC
         IWK(2*NFAC+I) = IWK(2*NFAC+I-1)*LFAC(NFAC-I+1)
   20 CONTINUE
      DO 340 I = 2, INTER
         IFAIL2 = 0
C
C        Compute number of terms for interaction using (COMBIN)
C
         NTERM = G02EAZ(I,NFAC,IFAIL2)
C
C        Set initial term
C
         DO 40 J = 1, I
            IWK(J) = 1
   40    CONTINUE
         DO 60 J = I + 1, NFAC
            IWK(J) = 0
   60    CONTINUE
C
C        Reject terms involving only main-plot terms
C
         DO 320 J = 1, NTERM
            IF (MPF.GT.0) THEN
               DO 80 K = MPF + 1, NFAC
                  IF (IWK(K).GT.0) GO TO 100
   80          CONTINUE
               GO TO 300
  100          CONTINUE
            END IF
C
            ITERM = ITERM + 1
C
C           ILEVEL = number of treatment combinations for term
C
            ILEVEL = 1
            NDF = 1
            DO 120 K = 1, NFAC
               IF (IWK(K).EQ.1) THEN
                  ILEVEL = ILEVEL*LFAC(K)
                  NDF = NDF*(LFAC(K)-1)
               END IF
  120       CONTINUE
            JJ = 0
            ICOMB = 0
C
C           Compute for first block
C
            DO 240 LL = 1, NTCPB
               DO 180 K = 1, ITC
                  ICOMB = ICOMB + 1
C
C                 Use G04CAU (SIMDO) to compute levels for each unit
C
                  CALL G04CAU(.TRUE.,.FALSE.,IWK(2*NFAC+1),NFAC,ICOMB,
     *                        IWK(NFAC+1),IFAIL2)
                  II = 1
                  KK = 1
C
C                 Compute combination index from levels
C
                  DO 140 L = NFAC, 1, -1
                     IF (IWK(L).EQ.1) THEN
                        II = II + KK*(IWK(NFAC+L)-1)
                        KK = KK*LFAC(L)
                     END IF
  140             CONTINUE
C
C                 Compute for replicate combinations
C
                  DO 160 L = 1, NREPL
                     JJ = JJ + 1
                     IFAC(JJ) = II
                     E(IE+II) = E(IE+II) + R(JJ)
                     TMEAN(IE+II) = TMEAN(IE+II) + Y(JJ)
  160             CONTINUE
  180          CONTINUE
C
C              Compute for replication due to main-plots
C
               DO 220 K = 2, MREP
                  DO 200 L = 1, KPLOT
                     JJ = JJ + 1
                     II = IFAC(JJ-KPLOT)
                     IFAC(JJ) = II
                     E(IE+II) = E(IE+II) + R(JJ)
                     TMEAN(IE+II) = TMEAN(IE+II) + Y(JJ)
  200             CONTINUE
  220          CONTINUE
  240       CONTINUE
C
C           Compute for remaining blocks
C
            DO 280 K = 2, NBLOCK
               DO 260 L = 1, KBLOCK
                  JJ = JJ + 1
                  II = IFAC(JJ-KBLOCK)
                  IFAC(JJ) = II
                  E(IE+II) = E(IE+II) + R(JJ)
                  TMEAN(IE+II) = TMEAN(IE+II) + Y(JJ)
  260          CONTINUE
  280       CONTINUE
C
C           Compute and remove effects
C
            IR = IR + 1
            CALL G04CAW(N,ILEVEL,E(IE+1),TMEAN(IE+1),SEMEAN(IR),ESS,R,
     *                  IFAC)
            TABLE(ITERM,1) = NDF
            TABLE(ITERM,2) = ESS
            NRDF = NRDF - NDF
            IE = IE + ILEVEL
            IMEAN(IR) = IE
  300       CONTINUE
C
C           G04CAV computes next factor combination
C
            CALL G04CAV(NFAC,I,IWK,IFAULT)
  320    CONTINUE
  340 CONTINUE
      RETURN
      END
