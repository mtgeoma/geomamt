      SUBROUTINE G08RBY(ZIN,ETA,VAPVEC,N,N1,NUNC,GAMMA,ICO,IRANK,IWA,
     *                  NIWA)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     ADJUSTS THE EXPECTED VALUES AND VARIANCE COVARIANCE MATRIX
C     OF A PARTICULAR FUNCTION OF THE ORDER STATISTICS FROM A
C     GENERALISED LOGISTIC DISTRIBUTION WHEN SOME OBSERVATIONS
C     ARE TIED.
C
C     PETTITT, A.N. APPROXIMATE METHODS USING RANKS FOR REGRESSION
C                   WITH CENSORED DATA.
C                   BIOMETRIKA, 70, PP 121-32.
C
C     ARGUMENTS :
C              ZIN - CONTAINS EXPECTED VALUES OF FUNCTION OF
C                    ORDER STATISTICS ON ENTRY AND ADJUSTED VALUES
C                    ON EXIT.
C              ETA - CONTAINS EXPECTED VALUES OF DERIVATIVE OF
C                    FUNCTION OF ORDER STATISTICS ON ENTRY AND
C                    ADJUSTED VALUES ON EXIT.
C           VAPVEC - CONTAINS VAR-COVAR MATRIX OF FUNCTION OF
C                    ORDER STATISTICS  ENTRY AND ADJUSTED MATRIX
C                    ON EXIT.
C                N - ON ENTRY SPECIFIES THE SAMPLE SIZE.
C               N1 - ON ENTRY SPECIFIES THE LENGTH OF A MATRIX
C                    REQUIRED TO STORE THE UPPER TRIANGLE OF A
C                    SQUARE SYMMETRIC MATRIX MATRIX OF ORDER N.
C                    (N1 .GE. N(N+1)/2).
C             NUNC - ON ENTRY SPECIFIES THE NUMBER OF UNCENSORED
C                    OBSERVATIONS IN THE SAMPLE.
C            GAMMA - ON ENTRY GAMMA SPECIFIES THE VALUE OF THE
C                    POWER PARAMETER IN THE GENERALISED LOGISTIC
C                    DISTRIBUTION.
C              ICO - ON ENTRY ICO SPECIFIES THE NUMBER OF SETS
C                    OF TIES.
C            IRANK - ON ENTRY IRANK CONTAINS THE RANKS OF THE
C                    OBSERVATIONS.
C              IWA - ON ENTRY, IWA CONTAINS THE FOLLOWING; AN INDEX
C                    TO THE ORIGINAL OBSERVATIONS, THE CENSORING OF THE
C                    ORIGINAL OBSERVATIONS(0 = UNCENSORED, 1 =
C                    CENSORED) THE NUMBER OF CENSORED OBSERVATIONS
C                    BETWEEN EACH PAIR OF CENSORED OBSERVATIONS AND
C                    THE POSITION AND LENGTH OF ANY SETS OF TIES.
C             NIWA - ON ENTRY, NIWA SPECIFIES THE LENGTH OF IWA.
C                    NIWA .GE. 4*N.
C
C
C     IF THERE ARE NO TIES, AMEND ZIN AND ETA AND RETURN
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  GAMMA
      INTEGER           ICO, N, N1, NIWA, NUNC
C     .. Array Arguments ..
      DOUBLE PRECISION  ETA(N), VAPVEC(N1), ZIN(N)
      INTEGER           IRANK(N), IWA(NIWA)
C     .. Local Scalars ..
      DOUBLE PRECISION  CNST, CRCT, EMOD, SUM1, SUM2, SUM3, SUM4, VMOD,
     *                  VMOD1, Z, ZMOD
      INTEGER           I, IFOUND, II, J, JC, JC1, JJ, JS, JS1, K, M,
     *                  N2, N3, NUN
C     .. External Functions ..
      INTEGER           G01DCU
      EXTERNAL          G01DCU
C     .. Executable Statements ..
      N2 = 2*N
      N3 = 3*N
      IF (ICO.EQ.0) THEN
         I = 0
   20    I = I + 1
         IF (IWA(N+IWA(I)).NE.0) GO TO 20
         DO 40 II = I, NUNC
            IF (GAMMA.LE.0.0001D0) THEN
               ZIN(IWA(I)) = ZIN(IWA(I)) - 1.0D0
            ELSE
               ZIN(IWA(I)) = (1-ZIN(IWA(I)))/GAMMA - ZIN(IWA(I))
               ETA(IWA(I)) = (1+GAMMA)*ETA(IWA(I))/GAMMA
            END IF
            I = I + IWA(N2+II) + 1
   40    CONTINUE
         RETURN
      END IF
C
C     IF ALL OBSERVATIONS ARE TIED RETURN
C
      IF (IWA(N3+2).EQ.NUNC) RETURN
C
      CNST = (1+1/GAMMA)*(1+1/GAMMA)
      IF (GAMMA.LE.0.0001D0) CNST = 1.0D0
      DO 440 JJ = 1, ICO
C
C        CHANGE MEANS AND VARIANCES OF EACH SET OF TIED
C        UNCENSORED OBSERVATIONS.
C
         JS = IWA(N3+2*(JJ-1)+1)
         JC = IWA(N3+2*(JJ-1)+2)
         SUM1 = 0.0D0
         SUM2 = 0.0D0
         SUM3 = 0.0D0
         SUM4 = 0.0D0
         DO 60 I = 1, JC
            Z = ZIN(IWA(JS+I-1))
            SUM1 = SUM1 + Z
            SUM2 = SUM2 + ETA(IWA(JS+I-1))
            M = G01DCU(IWA(JS+I-1),IWA(JS+I-1))
            SUM3 = SUM3 + VAPVEC(M)/CNST
            SUM4 = SUM4 + Z*Z
   60    CONTINUE
         CRCT = SUM4/JC - (SUM1/JC)*(SUM1/JC)
         IF (GAMMA.LE.0.0001D0) THEN
            ZMOD = SUM1/JC - 1.0D0
            EMOD = SUM2/JC
         ELSE
            ZMOD = (1-SUM1/JC)/GAMMA - SUM1/JC
            EMOD = (1+GAMMA)*(SUM2/JC)/GAMMA
         END IF
         VMOD = CNST*(SUM3/JC+CRCT)
         DO 80 I = 1, JC
            ZIN(IWA(JS+I-1)) = ZMOD
            ETA(IWA(JS+I-1)) = EMOD
            M = G01DCU(IWA(JS+I-1),IWA(JS+I-1))
            VAPVEC(M) = VMOD
   80    CONTINUE
C
C        CHANGE COVARIANCES OF TIED UNCENSORED OBSERVATIONS
C        IN THE SAME SET.
C
         SUM3 = 0.0D0
         DO 120 J = 1, JC
            DO 100 I = J + 1, JC
               M = G01DCU(IWA(JS+I-1),IWA(JS+J-1))
               SUM3 = SUM3 + VAPVEC(M)/CNST
  100       CONTINUE
  120    CONTINUE
         VMOD = CNST*(2*SUM3/(JC*(JC-1))-CRCT/(JC-1))
         DO 160 J = 1, JC
            DO 140 I = J + 1, JC
               M = G01DCU(IWA(JS+I-1),IWA(JS+J-1))
               VAPVEC(M) = VMOD
  140       CONTINUE
  160    CONTINUE
C
C        CHANGE COVARIANCES OF DIFFERENT SETS OF TIED UNCENSORED
C        OBSERVATIONS.
C
         DO 260 II = 1, JJ - 1
            JS1 = IWA(N3+2*(II-1)+1)
            JC1 = IWA(N3+2*(II-1)+2)
            SUM3 = 0.0D0
            DO 200 I = 1, JC1
               DO 180 J = 1, JC
                  M = G01DCU(IWA(JS1+I-1),IWA(JS+J-1))
                  SUM3 = SUM3 + VAPVEC(M)
  180          CONTINUE
  200       CONTINUE
            VMOD = SUM3/(JC*JC1)
            DO 240 I = 1, JC1
               DO 220 J = 1, JC
                  M = G01DCU(IWA(JS1+I-1),IWA(JS+J-1))
                  VAPVEC(M) = VMOD
  220          CONTINUE
  240       CONTINUE
  260    CONTINUE
C
C        CHANGE COVARIANCES OF TIED OBSERVATIONS WITH UNTIED - CENSORED
C        AND UNCENSORED.
C
         I = 0
         II = 1
         NUN = 0
  280    I = I + 1
         IFOUND = 0
         IF (I.GT.N) GO TO 440
         IF (II.GT.ICO) GO TO 300
         JS1 = IWA(N3+2*(II-1)+1)
         JC1 = IWA(N3+2*(II-1)+2)
  300    IF (I.EQ.JS1) THEN
            NUN = NUN + JC1
            I = I + JC1 - 1
            II = II + 1
            GO TO 280
         END IF
         IF (IWA(N+IWA(I)).EQ.0) THEN
            NUN = NUN + 1
            IFOUND = 1
         END IF
         IF (IFOUND.EQ.0) GO TO 340
C
C        HAVE NOW FOUND AN UNTIED UNCENSORED OBSERVATION
C
         IF (JJ.EQ.1) THEN
            IF (GAMMA.LE.0.0001D0) THEN
               ZIN(IWA(I)) = ZIN(IWA(I)) - 1.0D0
            ELSE
               ZIN(IWA(I)) = (1-ZIN(IWA(I)))/GAMMA - ZIN(IWA(I))
               ETA(IWA(I)) = (1+GAMMA)*ETA(IWA(I))/GAMMA
            END IF
         END IF
         VMOD = 0.0D0
         DO 320 J = 1, JC
            M = G01DCU(IWA(I),IWA(JS+J-1))
            VMOD = VMOD + VAPVEC(M)
  320    CONTINUE
  340    VMOD = VMOD/JC
         VMOD1 = 0.0D0
         IF (IWA(N2+NUN).EQ.0) GO TO 380
         DO 360 J = 1, JC
            M = G01DCU(IWA(I+IFOUND),IWA(JS+J-1))
            VMOD1 = VMOD1 + VAPVEC(M)
  360    CONTINUE
  380    VMOD1 = VMOD1/JC
         DO 420 J = 1, JC
            DO 400 K = 1, IFOUND + IWA(N2+NUN)
               M = G01DCU(IWA(I+K-1),IWA(JS+J-1))
               VAPVEC(M) = (1-IWA(N+IWA(I+K-1)))*VMOD + IWA(N+IWA(I+K-1)
     *                     )*VMOD1
  400       CONTINUE
  420    CONTINUE
         I = I + IFOUND + IWA(N2+NUN) - 1
         GO TO 280
  440 CONTINUE
      RETURN
      END
