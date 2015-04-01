      SUBROUTINE G02GCU(N,LINK,A,Y,FV,ETA,WT,NO,YMIN,Y0)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-937 (APR 1991).
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, Y0, YMIN
      INTEGER           N, NO
      CHARACTER         LINK
C     .. Array Arguments ..
      DOUBLE PRECISION  ETA(N), FV(N), WT(*), Y(N)
C     .. Local Scalars ..
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         LOG, SQRT
C     .. Executable Statements ..
      IF (N.EQ.NO) THEN
         IF (LINK.EQ.'I' .OR. LINK.EQ.'i') THEN
            DO 20 I = 1, N
               IF (Y(I).NE.0.0D0) THEN
                  FV(I) = Y(I)
                  ETA(I) = FV(I)
               ELSE
                  FV(I) = Y0
                  ETA(I) = Y0
               END IF
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
               IF (Y(I).GT.0.0D0) THEN
                  FV(I) = Y(I)
                  ETA(I) = SQRT(FV(I))
               ELSE IF (Y(I).EQ.0.0D0) THEN
                  FV(I) = Y0
                  ETA(I) = SQRT(Y0)
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
            DO 100 I = 1, N
               IF (Y(I).GT.0.0D0) THEN
                  FV(I) = Y(I)
                  ETA(I) = FV(I)**A
               ELSE
                  FV(I) = YMIN
                  ETA(I) = YMIN**A
               END IF
  100       CONTINUE
         END IF
      ELSE
         IF (LINK.EQ.'I' .OR. LINK.EQ.'i') THEN
            DO 120 I = 1, N
               IF (WT(I).GT.0.0D0) THEN
                  IF (Y(I).NE.0.0D0) THEN
                     FV(I) = Y(I)
                     ETA(I) = FV(I)
                  ELSE
                     FV(I) = Y0
                     ETA(I) = Y0
                  END IF
               ELSE
                  FV(I) = 0.0D0
                  ETA(I) = 0.0D0
               END IF
  120       CONTINUE
         ELSE IF (LINK.EQ.'L' .OR. LINK.EQ.'l') THEN
            DO 140 I = 1, N
               IF (WT(I).GT.0.0D0) THEN
                  IF (Y(I).GT.0.0D0) THEN
                     FV(I) = Y(I)
                     ETA(I) = LOG(FV(I))
                  ELSE
                     FV(I) = YMIN
                     ETA(I) = LOG(YMIN)
                  END IF
               ELSE
                  ETA(I) = 0.0D0
                  FV(I) = 0.0D0
               END IF
  140       CONTINUE
         ELSE IF (LINK.EQ.'S' .OR. LINK.EQ.'s') THEN
            DO 160 I = 1, N
               IF (WT(I).GT.0.0D0) THEN
                  IF (Y(I).GT.0.0D0) THEN
                     FV(I) = Y(I)
                     ETA(I) = SQRT(FV(I))
                  ELSE IF (Y(I).EQ.0.0D0) THEN
                     FV(I) = Y0
                     ETA(I) = SQRT(Y0)
                  ELSE
                     FV(I) = YMIN
                     ETA(I) = SQRT(YMIN)
                  END IF
               ELSE
                  ETA(I) = 0.0D0
                  FV(I) = 0.0D0
               END IF
  160       CONTINUE
         ELSE IF (LINK.EQ.'R' .OR. LINK.EQ.'r') THEN
            DO 180 I = 1, N
               IF (WT(I).GT.0.0D0) THEN
                  IF (Y(I).NE.0.0D0) THEN
                     FV(I) = Y(I)
                     ETA(I) = 1.0D0/FV(I)
                  ELSE
                     FV(I) = YMIN
                     ETA(I) = 1.0D0/YMIN
                  END IF
               ELSE
                  FV(I) = 0.0D0
                  ETA(I) = 0.0D0
               END IF
  180       CONTINUE
         ELSE IF (LINK.EQ.'E' .OR. LINK.EQ.'e') THEN
            DO 200 I = 1, N
               IF (WT(I).GT.0.0D0) THEN
                  IF (Y(I).GT.0.0D0) THEN
                     FV(I) = Y(I)
                     ETA(I) = FV(I)**A
                  ELSE
                     FV(I) = YMIN
                     ETA(I) = YMIN**A
                  END IF
               ELSE
                  FV(I) = 0.0D0
                  ETA(I) = 0.0D0
               END IF
  200       CONTINUE
         END IF
      END IF
      RETURN
      END
