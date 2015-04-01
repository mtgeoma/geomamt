      SUBROUTINE G02GCY(N,LINK,ETA,FV,T,A,WT,NO,IND)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     CALCULATES FITTED VALUES FROM LINEAR PREDICTORS
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A
      INTEGER           IND, N, NO
      CHARACTER         LINK
C     .. Array Arguments ..
      DOUBLE PRECISION  ETA(N), FV(N), T(*), WT(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AI, RAI, UFLO
      INTEGER           I, IAI
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      EXTERNAL          X02AMF
C     .. External Subroutines ..
      EXTERNAL          DCOPY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, EXP, LOG
C     .. Executable Statements ..
      IND = 1
      IF (LINK.EQ.'I' .OR. LINK.EQ.'i') THEN
         CALL DCOPY(N,ETA,1,FV,1)
      ELSE IF (LINK.EQ.'S' .OR. LINK.EQ.'s') THEN
         DO 20 I = 1, N
            FV(I) = ETA(I)*ETA(I)
   20    CONTINUE
      ELSE IF (N.EQ.NO) THEN
         IF (LINK.EQ.'L' .OR. LINK.EQ.'l') THEN
            UFLO = -LOG(X02AMF())
            DO 40 I = 1, N
               IF (ABS(ETA(I)).GT.UFLO) THEN
                  IND = 0
                  RETURN
               ELSE
                  FV(I) = EXP(ETA(I))
               END IF
   40       CONTINUE
         ELSE IF (LINK.EQ.'R' .OR. LINK.EQ.'r') THEN
            DO 60 I = 1, N
               IF (ETA(I).EQ.0.0D0) THEN
                  IND = 0
                  RETURN
               ELSE
                  FV(I) = 1.0D0/ETA(I)
               END IF
   60       CONTINUE
         ELSE IF (LINK.EQ.'E' .OR. LINK.EQ.'e') THEN
            AI = 1.0D0/A
            IAI = AI
            RAI = IAI
            IF (RAI.EQ.AI) THEN
               IF (IAI.GT.0) THEN
                  DO 80 I = 1, N
                     FV(I) = ETA(I)**IAI
   80             CONTINUE
               ELSE
                  DO 100 I = 1, N
                     IF (ETA(I).EQ.0.0D0) THEN
                        IND = 0
                        RETURN
                     ELSE
                        FV(I) = ETA(I)**IAI
                     END IF
  100             CONTINUE
               END IF
            ELSE IF (AI.GT.1.0D0) THEN
               DO 120 I = 1, N
                  IF (ETA(I).LT.0.0D0) THEN
                     IND = 0
                     RETURN
                  ELSE
                     FV(I) = ETA(I)**AI
                  END IF
  120          CONTINUE
            ELSE
               DO 140 I = 1, N
                  IF (ETA(I).LE.0.0D0) THEN
                     IND = 0
                     RETURN
                  ELSE
                     FV(I) = ETA(I)**AI
                  END IF
  140          CONTINUE
            END IF
         END IF
      ELSE
         IF (LINK.EQ.'L' .OR. LINK.EQ.'l') THEN
            UFLO = -LOG(X02AMF())
            DO 160 I = 1, N
               IF (WT(I).NE.0.0D0) THEN
                  IF (ABS(ETA(I)).GT.UFLO) THEN
                     IND = 0
                     RETURN
                  ELSE
                     FV(I) = EXP(ETA(I))
                  END IF
               ELSE
                  FV(I) = 0.0D0
               END IF
  160       CONTINUE
         ELSE IF (LINK.EQ.'R' .OR. LINK.EQ.'r') THEN
            DO 180 I = 1, N
               IF (WT(I).NE.0.0D0) THEN
                  IF (ETA(I).EQ.0.0D0) THEN
                     IND = 0
                     RETURN
                  ELSE
                     FV(I) = 1.0D0/ETA(I)
                  END IF
               ELSE
                  FV(I) = 0.0D0
               END IF
  180       CONTINUE
         ELSE IF (LINK.EQ.'E' .OR. LINK.EQ.'e') THEN
            AI = 1.0D0/A
            IAI = AI
            RAI = IAI
            IF (RAI.EQ.AI) THEN
               IF (IAI.GE.1) THEN
                  DO 200 I = 1, N
                     IF (WT(I).GT.0.0D0) THEN
                        FV(I) = ETA(I)**IAI
                     ELSE
                        FV(I) = 0.0D0
                     END IF
  200             CONTINUE
               ELSE
                  DO 220 I = 1, N
                     IF (WT(I).NE.0.0D0) THEN
                        IF (ETA(I).EQ.0.0D0) THEN
                           IND = 0
                           RETURN
                        ELSE
                           FV(I) = ETA(I)**IAI
                        END IF
                     ELSE
                        FV(I) = 0.0D0
                     END IF
  220             CONTINUE
               END IF
            ELSE IF (AI.GT.1.0D0) THEN
               DO 240 I = 1, N
                  IF (WT(I).NE.0.0D0) THEN
                     IF (ETA(I).LT.0.0D0) THEN
                        IND = 0
                        RETURN
                     ELSE
                        FV(I) = ETA(I)**AI
                     END IF
                  ELSE
                     FV(I) = 0.0D0
                  END IF
  240          CONTINUE
            ELSE
               DO 260 I = 1, N
                  IF (WT(I).NE.0.0D0) THEN
                     IF (ETA(I).LE.0.0D0) THEN
                        IND = 0
                        RETURN
                     ELSE
                        FV(I) = ETA(I)**AI
                     END IF
                  ELSE
                     FV(I) = 0.0D0
                  END IF
  260          CONTINUE
            END IF
         END IF
      END IF
      END
