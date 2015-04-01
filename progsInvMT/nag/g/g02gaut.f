      SUBROUTINE G02GAU(N,LINK,A,Y,FV,ETA,WT,NO,YMIN)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-920 (APR 1991).
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, YMIN
      INTEGER           N, NO
      CHARACTER*1       LINK
C     .. Array Arguments ..
      DOUBLE PRECISION  ETA(N), FV(N), WT(*), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  RIA
      INTEGER           I, IA
C     .. Intrinsic Functions ..
      INTRINSIC         LOG, SQRT
C     .. Executable Statements ..
      IF (N.EQ.NO) THEN
         IF (LINK.EQ.'I' .OR. LINK.EQ.'i') THEN
            DO 20 I = 1, N
               FV(I) = Y(I)
               ETA(I) = FV(I)
   20       CONTINUE
         ELSE IF (LINK.EQ.'L' .OR. LINK.EQ.'l') THEN
            DO 40 I = 1, N
               IF (Y(I).GT.0.0D0) THEN
                  FV(I) = Y(I)
                  ETA(I) = LOG(FV(I))
               ELSE
                  FV(I) = YMIN
                  ETA(I) = LOG(YMIN)
               END IF
   40       CONTINUE
         ELSE IF (LINK.EQ.'S' .OR. LINK.EQ.'s') THEN
            DO 60 I = 1, N
               IF (Y(I).GE.0.0D0) THEN
                  FV(I) = Y(I)
                  ETA(I) = SQRT(FV(I))
               ELSE
                  FV(I) = YMIN
                  ETA(I) = SQRT(YMIN)
               END IF
   60       CONTINUE
         ELSE IF (LINK.EQ.'R' .OR. LINK.EQ.'r') THEN
            DO 80 I = 1, N
               IF (Y(I).NE.0.0D0) THEN
                  FV(I) = Y(I)
                  ETA(I) = 1.0D0/FV(I)
               ELSE
                  FV(I) = YMIN
                  ETA(I) = 1.0D0/YMIN
               END IF
   80       CONTINUE
         ELSE IF (LINK.EQ.'E' .OR. LINK.EQ.'e') THEN
            IA = A
            RIA = IA
            IF (RIA.EQ.A) THEN
               IF (IA.EQ.1) THEN
                  DO 100 I = 1, N
                     ETA(I) = Y(I)
                     FV(I) = Y(I)
  100             CONTINUE
               ELSE
                  DO 120 I = 1, N
                     IF (Y(I).LE.0.0D0) THEN
                        ETA(I) = YMIN**IA
                        FV(I) = Y(I)
                     ELSE
                        ETA(I) = Y(I)**IA
                        FV(I) = Y(I)
                     END IF
  120             CONTINUE
               END IF
            ELSE
               DO 140 I = 1, N
                  IF (Y(I).LE.0.0D0) THEN
                     ETA(I) = YMIN**A
                     FV(I) = YMIN
                  ELSE
                     ETA(I) = Y(I)**A
                     FV(I) = Y(I)
                  END IF
  140          CONTINUE
            END IF
         END IF
      ELSE
         IF (LINK.EQ.'I' .OR. LINK.EQ.'i') THEN
            DO 160 I = 1, N
               IF (WT(I).GT.0.0D0) THEN
                  FV(I) = Y(I)
                  ETA(I) = FV(I)
               ELSE
                  FV(I) = 0.0D0
                  ETA(I) = 0.0D0
               END IF
  160       CONTINUE
         ELSE IF (LINK.EQ.'L' .OR. LINK.EQ.'l') THEN
            DO 180 I = 1, N
               IF (WT(I).GT.0.0D0) THEN
                  IF (Y(I).GT.0.0D0) THEN
                     FV(I) = Y(I)
                     ETA(I) = LOG(FV(I))
                  ELSE
                     FV(I) = YMIN
                     ETA(I) = LOG(YMIN)
                  END IF
               ELSE
                  FV(I) = 0.0D0
                  ETA(I) = 0.0D0
               END IF
  180       CONTINUE
         ELSE IF (LINK.EQ.'S' .OR. LINK.EQ.'s') THEN
            DO 200 I = 1, N
               IF (WT(I).GT.0.0D0) THEN
                  IF (Y(I).GE.0.0D0) THEN
                     FV(I) = Y(I)
                     ETA(I) = SQRT(FV(I))
                  ELSE
                     FV(I) = YMIN
                     ETA(I) = SQRT(YMIN)
                  END IF
               ELSE
                  ETA(I) = 0.0D0
                  FV(I) = 0.0D0
               END IF
  200       CONTINUE
         ELSE IF (LINK.EQ.'R' .OR. LINK.EQ.'r') THEN
            DO 220 I = 1, N
               IF (WT(I).GT.0.0D0) THEN
                  IF (Y(I).NE.0.0D0) THEN
                     FV(I) = Y(I)
                     ETA(I) = 1.0D0/FV(I)
                  ELSE
                     FV(I) = YMIN
                     ETA(I) = 1.0D0/YMIN
                  END IF
               ELSE
                  ETA(I) = 0.0D0
                  FV(I) = 0.0D0
               END IF
  220       CONTINUE
         ELSE IF (LINK.EQ.'E' .OR. LINK.EQ.'e') THEN
            IA = A
            RIA = IA
            IF (RIA.EQ.A) THEN
               IF (IA.EQ.1) THEN
                  DO 240 I = 1, N
                     IF (WT(I).GT.0.0D0) THEN
                        ETA(I) = Y(I)
                        FV(I) = Y(I)
                     ELSE
                        FV(I) = 0.0D0
                        ETA(I) = 0.0D0
                     END IF
  240             CONTINUE
               ELSE
                  DO 260 I = 1, N
                     IF (WT(I).EQ.0.0D0) THEN
                        ETA(I) = 0.0D0
                        FV(I) = 0.0D0
                     ELSE IF (Y(I).LE.0.0D0) THEN
                        ETA(I) = YMIN**IA
                        FV(I) = YMIN
                     ELSE
                        ETA(I) = Y(I)**IA
                        FV(I) = Y(I)
                     END IF
  260             CONTINUE
               END IF
            ELSE
               DO 280 I = 1, N
                  IF (WT(I).EQ.0.0D0) THEN
                     FV(I) = 0.0D0
                     ETA(I) = 0.0D0
                  ELSE IF (Y(I).LE.0.0D0) THEN
                     ETA(I) = YMIN**A
                     FV(I) = YMIN
                  ELSE
                     ETA(I) = Y(I)**A
                     FV(I) = Y(I)
                  END IF
  280          CONTINUE
            END IF
         END IF
      END IF
      END
