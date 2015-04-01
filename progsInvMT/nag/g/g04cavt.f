      SUBROUTINE G04CAV(N,M,IND,IERROR)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Generates subsets of M out of N in lexical order
C
C     IND is indicator array = 1 if in subset 0 otherwise
C
C     IND should be initialised to have only M items = 1.
C
C     .. Scalar Arguments ..
      INTEGER           IERROR, M, N
C     .. Array Arguments ..
      INTEGER           IND(N)
C     .. Local Scalars ..
      INTEGER           I, K, L
C     .. Executable Statements ..
      IF (N.LE.1) THEN
         IERROR = 1
      ELSE IF (M.GE.N) THEN
         IERROR = 1
      ELSE
         IERROR = 0
         K = 0
         DO 20 I = 1, N
            IF (IND(I).EQ.1) THEN
               K = K + 1
            ELSE IF (IND(I).NE.0) THEN
               IERROR = 2
               GO TO 180
            END IF
   20    CONTINUE
         IF (K.NE.M) THEN
            IERROR = 3
            GO TO 180
         END IF
C
C        Find last included item
C
         L = N
   40    CONTINUE
         IF (IND(L).EQ.0) THEN
            L = L - 1
            GO TO 40
         END IF
         IF (L.LT.N) THEN
C
C           If new last item swopped with following
C
            IND(L) = 0
            IND(L+1) = 1
         ELSE
C
C           if not find last item such that the following
C           can be included
C
            K = 1
            L = L - 1
   60       CONTINUE
            IF (IND(L).EQ.1) THEN
               L = L - 1
               K = K + 1
               GO TO 60
            END IF
            IF (K.LT.M) THEN
C
C              If sequence not complete
C
               L = L - 1
   80          CONTINUE
               IF (IND(L).EQ.0) THEN
                  L = L - 1
                  GO TO 80
               END IF
C
C              Swap last item with following
C
               IND(L) = 0
               L = L + 1
               IND(L) = 1
               IF (L+K.LT.N) THEN
C
C                 Include all following items
C
                  DO 100 I = L + 1, L + K
                     IND(I) = 1
  100             CONTINUE
                  DO 120 I = L + K + 1, N
                     IND(I) = 0
  120             CONTINUE
               END IF
            ELSE
C
C              Start again
C
               DO 140 I = 1, M
                  IND(I) = 1
  140          CONTINUE
               DO 160 I = M + 1, N
                  IND(I) = 0
  160          CONTINUE
            END IF
         END IF
      END IF
  180 CONTINUE
      RETURN
      END
