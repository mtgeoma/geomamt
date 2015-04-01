      SUBROUTINE G08RAV(Y,X,NX,IDIST,NS,NV,NMAX,N1,IP,NSUM,TOL,NPEST,
     *                  NPVAR,NWA,IRANK,ZIN,ETA,VAPVEC,PAREST,PARVAR,WA,
     *                  WB,IWA,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-640 (APR 1988).
C     MARK 14 REVISED. IER-794 (DEC 1989).
C
C     BLACK BOX ROUTINE FOR K-SAMPLE RANK ANALYSIS. MARGINAL LIKELIHOOD
C     FOR RANKS OF OBSERVATIONS IS APPROXIMATED USING ORDER STATISTICS
C     AND PARAMETER ESTIMATES OBTAINED FOR GENERAL LINEAR MODEL.
C     PETTITT A.N.P. INFERENCE FOR THE LINEAR MODEL USING A
C                    LIKELIHOOD BASED ON RANKS.
C                    JRSS B, 44, PP 234-243.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL
      INTEGER           IDIST, IFAIL, IP, N1, NMAX, NPEST, NPVAR, NS,
     *                  NSUM, NWA, NX
C     .. Array Arguments ..
      DOUBLE PRECISION  ETA(NMAX), PAREST(NPEST), PARVAR(NPVAR,IP),
     *                  VAPVEC(N1), WA(NWA,IP), WB(NMAX), X(NX,IP),
     *                  Y(NSUM), ZIN(NMAX)
      INTEGER           IRANK(NMAX), IWA(NMAX), NV(NS)
C     .. Local Scalars ..
      INTEGER           I, IERROR, II, J, JJ, L, MPOS, MPOS1, N11, NTIED
C     .. External Subroutines ..
      EXTERNAL          F06QHF, F06FBF, G08RAW, G08RAX, G08RAY, G08RAZ,
     *                  M01DAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      IERROR = 0
      CALL F06QHF('General',NPVAR,IP,0.0D0,0.0D0,PARVAR,NPVAR)
      CALL F06FBF(NPEST,0.0D0,PAREST,1)
C
C     FOR EACH SAMPLE SET UP SCORES AND CALCULATE CONTRIBUTION TO
C     LIKELIHOOD. CHECK WHETHER ANY OTHERS WITH SAME SAMPLE SIZE
C     TO AVOID REGENERATING SCORES.
C
      MPOS = 1
      NTIED = 0
      DO 220 II = 1, NS
         IF (NV(II).GT.0) THEN
            DO 40 I = 1, NV(II)
               L = MPOS - 1 + I
               WB(I) = Y(L)
               DO 20 J = 1, IP
                  WA(I,J) = X(L,J)
   20          CONTINUE
   40       CONTINUE
            N11 = NV(II)*(NV(II)+1)/2
            CALL M01DAF(WB,1,NV(II),'A',IRANK,IERROR)
            DO 60 I = 1, NV(II)
               ZIN(I) = WB(I)
   60       CONTINUE
            DO 80 I = 1, NV(II)
               WB(IRANK(I)) = ZIN(I)
   80       CONTINUE
            CALL G08RAZ(ZIN,ETA,VAPVEC,NV(II),N11,IDIST)
            CALL G08RAY(ZIN,ETA,VAPVEC,WB,NV(II),N11,TOL,IWA)
            IF (IWA(2).EQ.NV(II)) THEN
               NTIED = NTIED + 1
               GO TO 200
            END IF
            CALL G08RAX(ZIN,ETA,VAPVEC,NV(II),N11,IRANK,WA,NWA,PAREST,
     *                  IP,PARVAR,NPVAR)
         ELSE
            GO TO 200
         END IF
         IF (II.EQ.NS) GO TO 220
         MPOS1 = MPOS
         DO 180 JJ = II + 1, NS
            MPOS1 = MPOS1 + ABS(NV(JJ-1))
            IF (NV(JJ).EQ.NV(II)) THEN
               DO 100 I = 1, NV(JJ)
                  L = MPOS1 - 1 + I
                  WB(I) = Y(L)
                  WA(I,1) = ETA(I)
  100          CONTINUE
               CALL M01DAF(WB,1,NV(JJ),'A',IRANK,IERROR)
               DO 120 I = 1, NV(JJ)
                  ETA(I) = WB(I)
  120          CONTINUE
               DO 160 I = 1, NV(JJ)
                  WB(IRANK(I)) = ETA(I)
                  L = MPOS1 - 1 + I
                  DO 140 J = 1, IP
                     WA(I,J) = X(L,J)
  140             CONTINUE
  160          CONTINUE
               CALL G08RAY(ZIN,ETA,VAPVEC,WB,NV(JJ),N11,TOL,IWA)
               IF (IWA(2).EQ.NV(JJ)) THEN
                  NTIED = NTIED + 1
                  GO TO 200
               END IF
               CALL G08RAX(ZIN,ETA,VAPVEC,NV(JJ),N11,IRANK,WA,NWA,
     *                     PAREST,IP,PARVAR,NPVAR)
               NV(JJ) = -NV(JJ)
            END IF
  180    CONTINUE
  200    MPOS = MPOS + ABS(NV(II))
  220 CONTINUE
C
      DO 240 I = 1, NS
         NV(I) = ABS(NV(I))
  240 CONTINUE
      IF (NTIED.EQ.NS) THEN
         IERROR = 3
         GO TO 260
      END IF
C
C     MATRIX CAN NOW BE INVERTED TO PROVIDE ESTIMATES.
C
      CALL G08RAW(PAREST,NPEST,PARVAR,NPVAR,IP,WA,NWA,WB,IERROR)
C
  260 IFAIL = IERROR
      RETURN
      END
