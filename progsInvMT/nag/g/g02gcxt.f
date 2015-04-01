      SUBROUTINE G02GCX(N,LINK,ETA,T,DER,A,WT,NO)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     CALCULATES DERIVATIVE OF THE LINK FUNCTION
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A
      INTEGER           N, NO
      CHARACTER         LINK
C     .. Array Arguments ..
      DOUBLE PRECISION  DER(N), ETA(N), T(*), WT(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AI, RAI
      INTEGER           I, IAI
C     .. External Subroutines ..
      EXTERNAL          F06FBF
C     .. Intrinsic Functions ..
      INTRINSIC         EXP
C     .. Executable Statements ..
      IF (LINK.EQ.'I' .OR. LINK.EQ.'i') THEN
         CALL F06FBF(N,1.0D0,DER,1)
      ELSE IF (LINK.EQ.'S' .OR. LINK.EQ.'s') THEN
         DO 20 I = 1, N
            DER(I) = 2.0D0*ETA(I)
   20    CONTINUE
      ELSE IF (N.EQ.NO) THEN
         IF (LINK.EQ.'L' .OR. LINK.EQ.'l') THEN
            DO 40 I = 1, N
               DER(I) = EXP(ETA(I))
   40       CONTINUE
         ELSE IF (LINK.EQ.'R' .OR. LINK.EQ.'r') THEN
            DO 60 I = 1, N
               DER(I) = -1.0D0/(ETA(I)*ETA(I))
   60       CONTINUE
         ELSE IF (LINK.EQ.'E' .OR. LINK.EQ.'e') THEN
            AI = 1.0D0/A
            IAI = AI
            RAI = IAI
            IF (RAI.EQ.AI) THEN
               IAI = IAI - 1
               IF (IAI.EQ.0) THEN
                  CALL F06FBF(N,1.0D0,DER,1)
               ELSE
                  DO 80 I = 1, N
                     DER(I) = RAI*ETA(I)**IAI
   80             CONTINUE
               END IF
            ELSE
               RAI = AI - 1.0D0
               DO 100 I = 1, N
                  DER(I) = AI*ETA(I)**RAI
  100          CONTINUE
            END IF
         END IF
      ELSE
         IF (LINK.EQ.'L' .OR. LINK.EQ.'l') THEN
            DO 120 I = 1, N
               IF (WT(I).GT.0.0D0) THEN
                  DER(I) = EXP(ETA(I))
               ELSE
                  DER(I) = 0.0D0
               END IF
  120       CONTINUE
         ELSE IF (LINK.EQ.'R' .OR. LINK.EQ.'r') THEN
            DO 140 I = 1, N
               IF (WT(I).GT.0.0D0) THEN
                  DER(I) = -1.0D0/(ETA(I)*ETA(I))
               ELSE
                  DER(I) = 0.0D0
               END IF
  140       CONTINUE
         ELSE IF (LINK.EQ.'E' .OR. LINK.EQ.'e') THEN
            AI = 1.0D0/A
            IAI = AI
            RAI = IAI
            IF (RAI.EQ.AI) THEN
               IAI = IAI - 1
               IF (IAI.EQ.0) THEN
                  DO 160 I = 1, N
                     IF (WT(I).GT.0.0D0) THEN
                        DER(I) = 1.0D0
                     ELSE
                        DER(I) = 0.0D0
                     END IF
  160             CONTINUE
               ELSE
                  DO 180 I = 1, N
                     IF (WT(I).GT.0.0D0) THEN
                        DER(I) = RAI*ETA(I)**IAI
                     ELSE
                        DER(I) = 0.0D0
                     END IF
  180             CONTINUE
               END IF
            ELSE
               RAI = AI - 1.0D0
               DO 200 I = 1, N
                  IF (WT(I).GT.0.0D0) THEN
                     DER(I) = AI*ETA(I)**RAI
                  ELSE
                     DER(I) = 0.0D0
                  END IF
  200          CONTINUE
            END IF
         END IF
      END IF
      END
