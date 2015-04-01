      SUBROUTINE M01DBF(IV,M1,M2,ORDER,IRANK,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     M01DBF RANKS A VECTOR OF INTEGER NUMBERS IN ASCENDING
C     OR DESCENDING ORDER.
C
C     M01DBF USES A VARIANT OF LIST-MERGING, AS DESCRIBED
C     BY KNUTH. THE ROUTINE TAKES ADVANTAGE OF NATURAL
C     ORDERING IN THE DATA, AND USES A SIMPLE LIST INSERTION
C     IN A PREPARATORY PASS TO GENERATE ORDERED LISTS OF
C     LENGTH AT LEAST 10. THE RANKING IS STABLE: EQUAL ELEMENTS
C     PRESERVE THEIR ORDERING IN THE INPUT DATA.
C
C     THE MINIMUM LENGTH OF THE LISTS AT THE END OF THE
C     PREPARATORY PASS IS DEFINED BY THE VARIABLE MAXINS.
C
C     WRITTEN BY N.M.MACLAREN, UNIVERSITY OF CAMBRIDGE.
C     REVISED BY NAG CENTRAL OFFICE.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='M01DBF')
      INTEGER           MAXINS
      PARAMETER         (MAXINS=10)
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M1, M2
      CHARACTER*1       ORDER
C     .. Array Arguments ..
      INTEGER           IRANK(M2), IV(M2)
C     .. Local Scalars ..
      INTEGER           A, B, C, I, I1, I2, IERR, ILIST, J, J1, J2, K,
     *                  K1, K2, L, LIST1, LIST2, NLAST, NPREV, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
C
C       CHECK THE ARGUMENTS AND DEAL WITH THE TRIVIAL CASE.
C
      IF (M2.LT.1 .OR. M1.LT.1 .OR. M1.GT.M2) THEN
         IERR = 1
         WRITE (P01REC,FMT=99999) M1, M2
         NREC = 2
      ELSE IF (ORDER.NE.'A' .AND. ORDER.NE.'a' .AND. ORDER.NE.'D' .AND.
     *         ORDER.NE.'d') THEN
         IERR = 2
         WRITE (P01REC,FMT=99998) ORDER
         NREC = 1
      ELSE IF (M1.EQ.M2) THEN
         IRANK(M2) = M2
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
         K = M1
         I = K
         J = K
         L = K + MAXINS
         B = IV(K)
         C = B
         DO 40 K = M1 + 1, M2
C
C           DEAL WITH ADDITIONS AT EITHER END.
C
            A = IV(K)
            IF (A.GE.C) THEN
               IRANK(J) = K
               J = K
               C = A
            ELSE IF (A.LT.B) THEN
               IRANK(K) = I
               I = K
               B = A
            ELSE
C
C              DO AN ASCENDING LIST INSERTION.
C
               IF (K.LT.L) THEN
                  I2 = I
   20             I1 = I2
                  I2 = IRANK(I1)
                  IF (A.GE.IV(I2)) GO TO 20
                  IRANK(I1) = K
                  IRANK(K) = I2
               ELSE
C
C                 ADD THE CURRENT LIST ON TO THE OTHERS.
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
                  B = IV(K)
                  C = B
               END IF
            END IF
   40    CONTINUE
C
C        TIDY UP AT THE END.
C
         IRANK(J) = 0
         IF (ILIST.LT.0) THEN
            LIST1 = -I
            GO TO 280
         ELSE IF (ILIST.EQ.0) THEN
            LIST2 = -I
         ELSE
            IRANK(NPREV) = -I
         END IF
         IRANK(NLAST) = 0
C
C        AT THIS POINT:
C        LIST1 = - (INDEX OF LEAST ELEMENT IN THE FIRST LIST)
C        LIST2 = - (INDEX OF LEAST ELEMENT IN THE SECOND LIST)
C        FOR EACH K, IRANK(K) IS THE INDEX OF THE NEXT ELEMENT IN THE
C        CURRENT LIST, EXCEPT THAT, IF THERE IS NO SUCH ELEMENT,
C        IRANK(K) IS - (INDEX OF THE LEAST ELEMENT IN THE NEXT LIST
C        BUT 1)  OR 0 IF THERE IS NO SUCH LIST.
C
C        START MERGING LISTS BY PAIRS.
C
   60    ILIST = -1
         I = -LIST1
         J = -LIST2
   80    K = I
         IF (IV(I).GT.IV(J)) K = J
         IF (ILIST.LT.0) THEN
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
C        NEW VARIABLES I2, J2 AND K2 ARE USED INSTEAD OF I, J AND K
C        WITHIN THE INNERMOST BLOCK TO ENCOURAGE OPTIMISING COMPILERS TO
C        STORE THEM IN REGISTERS.
C         I2 POINTS TO THE CURRENT ELEMENT IN THE FIRST LIST
C         J2 POINTS TO THE CURRENT ELEMENT IN THE SECOND LIST
C         K2 POINTS TO THE CURRENT ELEMENT IN THE MERGED LIST
C
         I2 = I
         J2 = J
         IF (K.NE.I2) GO TO 140
  100    A = IV(J2)
         K2 = I2
  120    I2 = K2
         K2 = IRANK(I2)
         IF (K2.LE.0) GO TO 180
         IF (A.GE.IV(K2)) GO TO 120
         IRANK(I2) = J2
         I2 = K2
  140    A = IV(I2)
         K2 = J2
  160    J2 = K2
         K2 = IRANK(J2)
         IF (K2.LE.0) GO TO 200
         IF (A.GT.IV(K2)) GO TO 160
         IRANK(J2) = I2
         J2 = K2
         GO TO 100
C
C        ADD THE REMAINS OF ONE LIST TO THE OTHER.
C
  180    K = 1
         I1 = K2
         GO TO 220
  200    K = 2
         J1 = K2
  220    I = I2
         J = J2
         IF (K.EQ.1) THEN
C
C           FIRST LIST IS EXHAUSTED
C
            IRANK(I) = J
            I = -I1
            J1 = J
  240       J = J1
            J1 = IRANK(J)
            IF (J1.GT.0) GO TO 240
            L = J
            J = -J1
         ELSE
C
C           SECOND LIST IS EXHAUSTED
C
            IRANK(J) = I
            J = -J1
            I1 = I
  260       I = I1
            I1 = IRANK(I)
            IF (I1.GT.0) GO TO 260
            L = I
            I = -I1
         END IF
C
C        TIDY UP AND CARRY ON IF NOT FINISHED.
C
         IF ((I.NE.0) .AND. (J.NE.0)) GO TO 80
         IRANK(L) = 0
         K = I + J
         IF (ILIST.GT.0) THEN
            IRANK(NLAST) = -K
            GO TO 60
         ELSE IF (K.NE.0) THEN
            LIST2 = -K
            GO TO 60
         END IF
C
C        IF DESCENDING, REVERSE ALL POINTERS BETWEEN EQUALITY
C        BLOCKS.
C
  280    IF ((ORDER.EQ.'D') .OR. (ORDER.EQ.'d')) THEN
            I = 0
            J = -LIST1
  300       K = J
            K1 = K
            A = IV(K)
  320       K = K1
            K1 = IRANK(K)
            IF (K1.NE.0) THEN
               IF (A.EQ.IV(K1)) GO TO 320
            END IF
            IRANK(K) = I
            I = J
            J = K1
            IF (J.NE.0) GO TO 300
            LIST1 = -I
         END IF
C
C        CONVERT THE LIST FORM TO RANKS AND RETURN.
C
         K = M1
         I = -LIST1
  340    I1 = IRANK(I)
         IRANK(I) = K
         K = K + 1
         I = I1
         IF (I.GT.0) GO TO 340
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
     *  'lues is illegal',/'    M1 =',I16,'  M2 =',I16)
99998 FORMAT (' ** On entry, ORDER has an illegal value: ORDER = ',A1)
      END
