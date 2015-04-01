      SUBROUTINE M01DJF(RM,LDM,M1,M2,N1,N2,ORDER,IRANK,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     M01DJF RANKS THE COLUMNS OF A MATRIX OF REAL NUMBERS IN
C     ASCENDING OR DESCENDING ORDER.
C
C     M01DJF RANKS COLUMNS N1 TO N2 OF A MATRIX, USING THE DATA
C     IN ROWS M1 TO M2 OF THOSE COLUMNS. THE ORDERING IS
C     DETERMINED BY FIRST RANKING THE DATA IN ROW M1, THEN
C     RANKING ANY TIED COLUMNS ACCORDING TO THE DATA IN ROW
C     M1+1, AND SO ON UP TO ROW M2.
C
C     M01DJF USES A VARIANT OF LIST-MERGING, AS DECRIBED BY
C     KNUTH. THE ROUTINE TAKES ADVANTAGE OF ANY NATURAL ORDERING
C     IN THE DATA, AND USES A SIMPLE LIST INSERTION IN A
C     PREPARATORY PASS TO GENERATE ORDERED LISTS OF LENGTH AT
C     LEAST 10. THE RANKING IS STABLE: EQUAL COLUMNS PRESERVE
C     THEIR ORDERING IN THE INPUT DATA.
C
C     THE MINIMUM LENGTH OF THE LISTS AT THE END OF THE
C     PREPARATORY PASS IS DEFINED BY THE VARIABLE MAXINS.
C
C     WRITTEN BY N.M.MACLAREN, UNIVERSITY OF CAMBRIDGE.
C     REVISED BY NAG CENTRAL OFFICE.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='M01DJF')
      INTEGER           MAXINS
      PARAMETER         (MAXINS=10)
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDM, M1, M2, N1, N2
      CHARACTER*1       ORDER
C     .. Array Arguments ..
      DOUBLE PRECISION  RM(LDM,N2)
      INTEGER           IRANK(N2)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, C
      INTEGER           I, I1, I2, IERR, IH, ILIST, J, J1, K, K1, KTEMP,
     *                  L, LIST1, LIST2, NLAST, NPREV, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(3)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
C
C       CHECK THE ARGUMENTS AND DEAL WITH THE TRIVIAL CASE.
C
      IF (M2.LT.1 .OR. N2.LT.1 .OR. M1.LT.1 .OR. M1.GT.M2 .OR. N1.LT.
     *    1 .OR. N1.GT.N2 .OR. LDM.LT.M2) THEN
         IERR = 1
         WRITE (P01REC,FMT=99999) M1, M2, LDM, N1, N2
         NREC = 3
      ELSE IF (ORDER.NE.'A' .AND. ORDER.NE.'a' .AND. ORDER.NE.'D' .AND.
     *         ORDER.NE.'d') THEN
         IERR = 2
         WRITE (P01REC,FMT=99998) ORDER
         NREC = 1
      ELSE IF (N1.EQ.N2) THEN
         IRANK(N2) = N2
         IERR = 0
      ELSE
         IERR = 0
C
C        INITIALISE, USING NATURAL RUNS IN BOTH DIRECTIONS AND
C        STRAIGHT LIST INSERTION FOR SMALL LISTS.
C
C        I  POINTS TO THE SMALLEST ELEMENT IN THE CURRENT LIST
C        J  POINTS TO THE LARGEST  ELEMENT IN THE CURRENT LIST
C        B  IS THE VALUE OF THE SMALLEST ELEMENT IN CURRENT LIST
C        C  IS THE VALUE OF THE LARGEST  ELEMENT IN CURRENT LIST
C
         ILIST = -1
         K = N1
         I = K
         J = K
         L = K + MAXINS
         B = RM(M1,K)
         C = B
         DO 200 K = N1 + 1, N2
C
C           DEAL WITH ADDITIONS AT EITHER END.
C
            A = RM(M1,K)
            IF (A.GE.C) THEN
               IF (A.EQ.C) THEN
                  DO 20 IH = M1 + 1, M2
                     IF (RM(IH,K).GT.RM(IH,J)) GO TO 40
                     IF (RM(IH,K).LT.RM(IH,J)) GO TO 60
   20             CONTINUE
               END IF
   40          IRANK(J) = K
               J = K
               C = A
               GO TO 200
            END IF
   60       IF (A.LE.B) THEN
               IF (A.EQ.B) THEN
                  DO 80 IH = M1 + 1, M2
                     IF (RM(IH,K).LT.RM(IH,I)) GO TO 100
                     IF (RM(IH,K).GT.RM(IH,I)) GO TO 120
   80             CONTINUE
                  GO TO 120
               END IF
  100          IRANK(K) = I
               I = K
               B = A
               GO TO 200
            END IF
C
C           DO AN ASCENDING LIST INSERTION.
C
  120       IF (K.LT.L) THEN
               I2 = I
  140          I1 = I2
               I2 = IRANK(I1)
               IF (A.GT.RM(M1,I2)) GO TO 140
               IF (A.EQ.RM(M1,I2)) THEN
                  DO 160 IH = M1 + 1, M2
                     IF (RM(IH,K).LT.RM(IH,I2)) GO TO 180
                     IF (RM(IH,K).GT.RM(IH,I2)) GO TO 140
  160             CONTINUE
                  GO TO 140
               END IF
  180          CONTINUE
               IRANK(I1) = K
               IRANK(K) = I2
            ELSE
C
C              ADD THE CURRENT LIST ON TO THE OTHERS.
C
               IF (ILIST.LT.0) THEN
                  LIST1 = -I
                  ILIST = 0
               ELSE IF (ILIST.EQ.0) THEN
                  LIST2 = -I
                  ILIST = 1
                  NPREV = NLAST
               ELSE
                  IRANK(NPREV) = -I
                  NPREV = NLAST
               END IF
C
               NLAST = J
               I = K
               J = K
               L = K + MAXINS
               B = RM(M1,K)
               C = B
            END IF
  200    CONTINUE
C
C        TIDY UP AT THE END.
C
         IRANK(J) = 0
         IF (ILIST.LT.0) THEN
            LIST1 = -I
            GO TO 580
         ELSE IF (ILIST.EQ.0) THEN
            LIST2 = -I
         ELSE
            IRANK(NPREV) = -I
         END IF
         IRANK(NLAST) = 0
C
C        AT THIS POINT:
C        LIST1  = -(INDEX OF LEAST ELEMENT IN FIRST LIST)
C        LIST2  = -(INDEX OF LEAST ELEMENT IN SECOND LIST)
C        FOR EACH K, IRANK(K) = INDEX OF NEXT ELEMENT IN CURRENT LIST,
C        EXCEPT THAT, IF THERE IS NO SUCH ELEMENT, IRANK(K) =
C        -(INDEX OF LEAST ELEMENT IN NEXT LIST BUT 1)  OR 0 IF THERE IS
C        NO SUCH LIST.
C
C
C        START MERGING LISTS BY PAIRS.
C
  220    ILIST = -1
         I = -LIST1
         J = -LIST2
  240    K = I
         IF (RM(M1,I).GT.RM(M1,J)) K = J
         IF (RM(M1,I).EQ.RM(M1,J)) THEN
            DO 260 IH = M1 + 1, M2
               IF (RM(IH,I).GT.RM(IH,J)) GO TO 280
               IF (RM(IH,I).LT.RM(IH,J)) GO TO 300
  260       CONTINUE
            GO TO 300
  280       K = J
         END IF
  300    IF (ILIST.LT.0) THEN
            LIST1 = -K
            ILIST = 0
         ELSE IF (ILIST.EQ.0) THEN
            LIST2 = -K
            ILIST = 1
            NLAST = L
         ELSE
            IRANK(NLAST) = -K
            NLAST = L
         END IF
C
C        MOVE ALONG THE LISTS UNTIL ONE FINISHES.
C
C        I  POINTS TO THE CURRENT ELEMENT IN THE FIRST LIST
C        J  POINTS TO THE CURRENT ELEMENT IN THE SECOND LIST
C        K  POINTS TO THE CURRENT ELEMENT IN THE MERGED LIST
C
         IF (K.NE.I) GO TO 400
  320    A = RM(M1,J)
         K = I
  340    I = K
         K = IRANK(I)
         IF (K.LE.0) GO TO 480
         IF (A.GT.RM(M1,K)) GO TO 340
         IF (A.EQ.RM(M1,K)) THEN
            DO 360 IH = M1 + 1, M2
               IF (RM(IH,J).LT.RM(IH,K)) GO TO 380
               IF (RM(IH,J).GT.RM(IH,K)) GO TO 340
  360       CONTINUE
            GO TO 340
         END IF
  380    IRANK(I) = J
         I = K
  400    A = RM(M1,I)
         K = J
  420    J = K
         K = IRANK(J)
         IF (K.LE.0) GO TO 500
         IF (A.GT.RM(M1,K)) GO TO 420
         IF (A.EQ.RM(M1,K)) THEN
            DO 440 IH = M1 + 1, M2
               IF (RM(IH,I).LT.RM(IH,K)) GO TO 460
               IF (RM(IH,I).GT.RM(IH,K)) GO TO 420
  440       CONTINUE
         END IF
  460    IRANK(J) = I
         J = K
         GO TO 320
C
C        ADD THE REMAINS OF ONE LIST TO THE OTHER.
C
  480    KTEMP = 1
         GO TO 520
  500    KTEMP = 2
  520    IF (KTEMP.EQ.1) THEN
C
C           FIRST LIST IS EXHAUSTED
C
            IRANK(I) = J
            I = -K
            J1 = J
  540       J = J1
            J1 = IRANK(J)
            IF (J1.GT.0) GO TO 540
            L = J
            J = -J1
         ELSE
C
C           SECOND LIST IS EXHAUSTED
C
            IRANK(J) = I
            J = -K
            I1 = I
  560       I = I1
            I1 = IRANK(I)
            IF (I1.GT.0) GO TO 560
            L = I
            I = -I1
         END IF
C
C        TIDY UP AND CARRY ON IF NOT FINISHED.
C
         IF ((I.NE.0) .AND. (J.NE.0)) GO TO 240
         IRANK(L) = 0
         K = I + J
         IF (ILIST.GT.0) THEN
            IRANK(NLAST) = -K
            GO TO 220
         ELSE IF (K.NE.0) THEN
            LIST2 = -K
            GO TO 220
         END IF
C
C        IF DESCENDING, REVERSE ALL POINTERS BETWEEN EQUALITY
C        BLOCKS.
C
  580    IF (ORDER.EQ.'D' .OR. ORDER.EQ.'d') THEN
            I = 0
            J = -LIST1
  600       K = J
            K1 = K
            A = RM(M1,K)
            KTEMP = K
  620       K = K1
            K1 = IRANK(K)
            IF (K1.NE.0) THEN
               IF (A.EQ.RM(M1,K1)) THEN
                  DO 640 IH = M1 + 1, M2
                     IF (RM(IH,KTEMP).NE.RM(IH,K1)) GO TO 660
  640             CONTINUE
                  GO TO 620
               END IF
            END IF
  660       IRANK(K) = I
            I = J
            J = K1
            IF (J.NE.0) GO TO 600
            LIST1 = -I
         END IF
C
C        CONVERT THE LIST FORM TO RANKS AND RETURN.
C
         K = N1
         I = -LIST1
  680    I1 = IRANK(I)
         IRANK(I) = K
         K = K + 1
         I = I1
         IF (I.GT.0) GO TO 680
C
      END IF
C
      IF (IERR.NE.0) THEN
         IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,P01REC)
      ELSE
         IFAIL = 0
      END IF
      RETURN
C
99999 FORMAT (' ** On entry, one or more of the following parameter va',
     *  'lues is illegal',/'    M1 =',I16,'  M2 =',I16,'  LDM =',I16,
     *  /'    N1 =',I16,'  N2 =',I16)
99998 FORMAT (' ** On entry, ORDER has an illegal value: ORDER = ',A1)
      END
