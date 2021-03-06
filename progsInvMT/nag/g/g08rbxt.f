      SUBROUTINE G08RBX(Y,ICEN,X,NX,GAMMA,NS,NV,NMAX,N1,IP,NSUM,TOL,
     *                  NPEST,NPVAR,NWA,NIWA,IRANK,ZIN,ETA,VAPVEC,
     *                  PAREST,PARVAR,WA,WB,IWA,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-641 (APR 1988).
C     MARK 14 REVISED. IER-795 (DEC 1989).
C
C     BLACK-BOX ROUTINE FOR K-SAMPLE RANK ANALYSIS WHEN SOME
C     OBSERVATIONS ARE RIGHT CENSORED. MARGINAL LIKELIHOOD FOR
C     RANKS OF OBSERVATIONS ASSUMING GENERALISED LOGISTIC
C     DISTRIBUTION IS APPROXIMATED USING ORDER STATISTICS AND
C     PARAMETER ESTIMATES OBTAINED FOR LINEAR MODEL.
C
C     PETTITT A.N.  APPROXIMATE METHODS USING RANKS FOR REGRESSION
C                   WITH CENSORED DATA.
C                   BIOMETRIKA, 70, PP 121-32.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  GAMMA, TOL
      INTEGER           IFAIL, IP, N1, NIWA, NMAX, NPEST, NPVAR, NS,
     *                  NSUM, NWA, NX
C     .. Array Arguments ..
      DOUBLE PRECISION  ETA(NMAX), PAREST(NPEST), PARVAR(NPVAR,IP),
     *                  VAPVEC(N1), WA(NWA,IP), WB(NMAX), X(NX,IP),
     *                  Y(NSUM), ZIN(NMAX)
      INTEGER           ICEN(NSUM), IRANK(NMAX), IWA(NIWA), NV(NS)
C     .. Local Scalars ..
      INTEGER           I, ICO, IERROR, II, JJ, L, N11, NSTART, NTIED,
     *                  NUNC
C     .. External Subroutines ..
      EXTERNAL          F06QHF, F06FBF, G08RAW, G08RAX, G08RBY, G08RBZ,
     *                  M01DAF
C     .. Executable Statements ..
      IERROR = 0
      CALL F06QHF('General',NPVAR,IP,0.0D0,0.0D0,PARVAR,NPVAR)
      CALL F06FBF(NPEST,0.0D0,PAREST,1)
C
C     FOR EACH SAMPLE CALCULATE CONTRIBUTION TO LIKELIHOOD AND
C     ADD TO TOTAL
C
      NSTART = 1
      NTIED = 0
      DO 140 I = 1, NS
         NUNC = 0
         DO 40 II = 1, NV(I)
            L = NSTART - 1 + II
            WB(II) = Y(L)
            DO 20 JJ = 1, IP
               WA(II,JJ) = X(L,JJ)
   20       CONTINUE
            IF (ICEN(L).EQ.0) NUNC = NUNC + 1
   40    CONTINUE
C
C        RANK EACH SET OF OBSERVATIONS
C
         CALL M01DAF(WB,1,NV(I),'A',IRANK,IERROR)
         DO 60 II = 1, NV(I)
            ZIN(II) = WB(II)
   60    CONTINUE
         DO 80 II = 1, NV(I)
            WB(IRANK(II)) = ZIN(II)
            IWA(IRANK(II)) = II
   80    CONTINUE
         DO 100 II = 1, NV(I)
            L = NSTART - 1 + II
            IWA(NV(I)+II) = ICEN(L)
  100    CONTINUE
C
C        CALCULATE SCORES FOR EACH SAMPLE
C
         N11 = NV(I)*(NV(I)+1)/2
         CALL G08RBZ(WB,NV(I),IRANK,NUNC,GAMMA,TOL,ZIN,ETA,VAPVEC,N11,
     *               ICO,IWA,NIWA)
C
C        AMEND SCORES IN CASE OF TIES
C
         CALL G08RBY(ZIN,ETA,VAPVEC,NV(I),N11,NUNC,GAMMA,ICO,IRANK,IWA,
     *               NIWA)
         IF (IWA(3*NV(I)+2).EQ.NUNC) NTIED = NTIED + 1
C
C        CALCULATE CONTRIBUTION TO LIKELIHOOD
C
         DO 120 II = 1, NV(I)
            IWA(II) = II
  120    CONTINUE
         CALL G08RAX(ZIN,ETA,VAPVEC,NV(I),N11,IWA,WA,NWA,PAREST,IP,
     *               PARVAR,NPVAR)
         NSTART = NSTART + NV(I)
  140 CONTINUE
      IF (NTIED.EQ.NS) THEN
         IERROR = 3
         GO TO 160
      END IF
C
C     CALCULATE PARAMETER ESTIMATES
C
      CALL G08RAW(PAREST,NPEST,PARVAR,NPVAR,IP,WA,NWA,WB,IERROR)
      IF (IERROR.NE.0) IERROR = 4
C
  160 IFAIL = IERROR
      RETURN
      END
