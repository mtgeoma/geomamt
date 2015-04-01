      SUBROUTINE Y90RNF(N,A,IA,NNVSTR,NVSTR,VSTR,STACK,PATH,VERTEX,
     *                  NUMBER,LOWLNK,NEDGE,LOEDGE)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ==========================================
C         *  Y90RNF :  Generate strong components  *
C         ==========================================
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
C     .. Scalar Arguments ..
      INTEGER           IA, N, NNVSTR
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*)
      INTEGER           LOEDGE(*), LOWLNK(*), NEDGE(*), NUMBER(*),
     *                  NVSTR(*), PATH(*), STACK(*), VERTEX(*), VSTR(*)
C     .. Local Scalars ..
      INTEGER           I, IFAIL, IPATH, ISTACK, IVCOMP, J, K, K1, K2,
     *                  L, L1, NN
      LOGICAL           LOOP2, LOOP3
C     .. External Subroutines ..
      EXTERNAL          M01DBF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     1. Generate strong components
C
C-----------------------------------------------------------------------
C
C     Initialize
C
      DO 40 I = 1, N
         VERTEX(I) = 0
         VSTR(I) = 0
         NN = 0
         DO 20 J = 1, N
            IF ((I.NE.J) .AND. (A(I,J).NE.ZERO)) THEN
               NN = NN + 1
            END IF
   20    CONTINUE
         NEDGE(I) = NN
         LOEDGE(I) = 1
   40 CONTINUE
      I = 0
      ISTACK = 0
      IPATH = 0
      IVCOMP = 0
      NNVSTR = 0
C
C     Main loop
C
      DO 280 L1 = 1, N
C
C     Select starting vertex
C
         K = 0
         DO 60 J = 1, N
            IF (VERTEX(J).LE.0) THEN
               K = J
               GO TO 80
            END IF
   60    CONTINUE
   80    CONTINUE
         LOOP2 = .TRUE.
C
C        New starting vertex was found
C
         IF (K.NE.0) THEN
C
  100       CONTINUE
            IF (LOOP2) THEN
C
C              Step 2 :  visit a vertex
C
               VERTEX(K) = 1
               I = I + 1
               NUMBER(K) = I
               LOWLNK(K) = I
               IPATH = IPATH + 1
               ISTACK = ISTACK + 1
               PATH(IPATH) = K
               STACK(ISTACK) = K
               LOOP3 = .TRUE.
C
  120          CONTINUE
               IF (LOOP3) THEN
C
C                    Step 3 :  explore an edge
C
                  K = PATH(IPATH)
                  IF (NEDGE(K).GE.1) THEN
                     NEDGE(K) = NEDGE(K) - 1
                     L = -1
                     DO 140 J = LOEDGE(K), N
                        IF ((J.NE.K) .AND. (A(K,J).NE.ZERO)) THEN
                           L = J
                           GO TO 160
                        END IF
  140                CONTINUE
  160                CONTINUE
                     LOEDGE(K) = L + 1
C
C                          Step 3.a
C
                     IF (VERTEX(L).LE.0) THEN
                        VERTEX(L) = 1
                        K = L
                        LOOP3 = .FALSE.
C
C                          Step 3.b and 4
C
                     ELSE
                        IF (NUMBER(L).LT.NUMBER(K)) THEN
                           DO 180 J = 1, ISTACK
                              IF (STACK(J).EQ.L) THEN
                                 LOWLNK(K) = MIN(LOWLNK(K),LOWLNK(L))
                                 GO TO 200
                              END IF
  180                      CONTINUE
  200                      CONTINUE
                        END IF
                        LOOP3 = .TRUE.
                     END IF
                  ELSE
C
C                          Step 3.c and 5
C
                     IF (LOWLNK(K).LT.NUMBER(K)) THEN
                        IPATH = IPATH - 1
                        L = PATH(IPATH)
                        LOWLNK(L) = MIN(LOWLNK(K),LOWLNK(L))
                        K = L
                        LOOP3 = .TRUE.
C
C                          Step 3.d and 6
C
                     ELSE
                        DO 220 J = 1, ISTACK
                           IF (STACK(J).EQ.K) THEN
                              L = J
                              GO TO 240
                           END IF
  220                   CONTINUE
  240                   CONTINUE
                        NNVSTR = NNVSTR + 1
                        DO 260 J = L, ISTACK
                           IVCOMP = IVCOMP + 1
                           VSTR(IVCOMP) = STACK(J)
  260                   CONTINUE
                        NVSTR(NNVSTR) = ISTACK - L + 1
                        ISTACK = L - 1
                        IPATH = IPATH - 1
                        IF (IPATH.GE.1) THEN
                           LOOP3 = .TRUE.
                        ELSE
                           LOOP2 = .FALSE.
                           LOOP3 = .FALSE.
                        END IF
                     END IF
                  END IF
               END IF
               IF (LOOP3) GO TO 120
C
            END IF
            IF (LOOP2) GO TO 100
C
         ELSE
            GO TO 300
         END IF
  280 CONTINUE
  300 CONTINUE
C-----------------------------------------------------------------------
C
C     2. Rank strong components in ascending order
C
C-----------------------------------------------------------------------
      IFAIL = 1
      CALL M01DBF(NVSTR,1,NNVSTR,'Ascending',PATH(1),IFAIL)
      DO 320 I = 1, NNVSTR
         STACK(PATH(I)) = I
  320 CONTINUE
C
      K2 = NVSTR(1)
      NVSTR(1) = 1
      DO 340 I = 1, NNVSTR
         K1 = K2
         K2 = NVSTR(I+1)
         NVSTR(I+1) = NVSTR(I) + K1
  340 CONTINUE
C-----------------------------------------------------------------------
C
C     End of Y90RNF
C
C-----------------------------------------------------------------------
      RETURN
      END
