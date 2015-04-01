      SUBROUTINE M01DZF(COMPAR,M1,M2,IRANK,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     M01DZF RANKS ARBITRARY DATA ACCORDING TO A USER-SUPPLIED
C     COMPARISON ROUTINE.
C
C     M01DZF IS A GENERAL-PURPOSE ROUTINE FOR RANKING ARBITRARY
C     DATA. M01DZF DOES NOT ACCESS THE DATA DIRECTLY; INSTEAD IT
C     CALLS A USER-SUPPLIED ROUTINE COMPAR TO DETERMINE THE RELATIVE
C     ORDERING OF ANY TWO DATA ITEMS. THE DATA ITEMS ARE IDENTIFIED
C     SIMPLY BY AN INTEGER IN THE RANGE M1 TO M2.
C
C     M01DZF USES A VARIANT OF LIST-MERGING, AS DESCRIBED BY KNUTH.
C     THE ROUTINE TAKES ADVANTAGE OF ANY NATURAL ORDERING IN THE
C     DATA, AND USES A SIMPLE LIST INSERTION IN A PREPARATORY PASS
C     TO GENERATE ORDERED LISTS OF LENGTH AT LEAST 10.
C
C     THE MINIMUM LENGTH OF THE LISTS AT THE END OF THE
C     PREPARATORY PASS IS DEFINED BY THE VARIABLE MAXINS.
C
C     WRITTEN BY NAG CENTRAL OFFICE.
C
C     .. Parameters ..
      INTEGER           MAXINS
      PARAMETER         (MAXINS=10)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='M01DZF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M1, M2
C     .. Array Arguments ..
      INTEGER           IRANK(M2)
C     .. Function Arguments ..
      LOGICAL           COMPAR
      EXTERNAL          COMPAR
C     .. Local Scalars ..
      INTEGER           I, I1, I2, IERR, ILIST, J, J1, K, KTEMP, L,
     *                  LIST1, LIST2, NLAST, NPREV
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
C
      ELSE IF (M1.EQ.M2) THEN
         IRANK(M1) = M1
         IERR = 0
      ELSE
         IERR = 0
C
C        INITIALISE, USING NATURAL RUNS IN BOTH DIRECTIONS AND
C        STRAIGHT LIST INSERTION FOR SMALL LISTS.
C
C        I  POINTS TO THE SMALLEST ELEMENT IN THE CURRENT LIST
C        J  POINTS TO THE LARGEST  ELEMENT IN THE CURRENT LIST
C
         ILIST = -1
         I = M1
         J = I
         L = I + MAXINS
         DO 40 K = M1 + 1, M2
C
C           DEAL WITH ADDITIONS AT EITHER END.
C
            IF ( .NOT. COMPAR(J,K)) THEN
               IRANK(J) = K
               J = K
            ELSE IF (COMPAR(I,K)) THEN
               IRANK(K) = I
               I = K
            ELSE
C
C              DO AN ASCENDING LIST INSERTION.
C
               IF (K.LT.L) THEN
                  I2 = I
   20             I1 = I2
                  I2 = IRANK(I1)
                  IF (COMPAR(K,I2)) GO TO 20
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
         IF (COMPAR(I,J)) K = J
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
C         I POINTS TO THE CURRENT ELEMENT IN THE FIRST LIST
C         J POINTS TO THE CURRENT ELEMENT IN THE SECOND LIST
C         K POINTS TO THE CURRENT ELEMENT IN THE MERGED LIST
C
         IF (K.NE.I) GO TO 140
  100    K = I
  120    I = K
         K = IRANK(I)
         IF (K.LE.0) GO TO 180
         IF (COMPAR(J,K)) GO TO 120
         IRANK(I) = J
         I = K
  140    K = J
  160    J = K
         K = IRANK(J)
         IF (K.LE.0) GO TO 200
         IF (COMPAR(I,K)) GO TO 160
         IRANK(J) = I
         J = K
         GO TO 100
C
C        ADD THE REMAINS OF ONE LIST TO THE OTHER.
C
  180    KTEMP = 1
         GO TO 220
  200    KTEMP = 2
  220    IF (KTEMP.EQ.1) THEN
C
C           FIRST LIST IS EXHAUSTED
C
            IRANK(I) = J
            I = -K
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
            J = -K
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
C        CONVERT THE LIST FORM TO RANKS AND RETURN.
C
  280    K = M1
         I = -LIST1
  300    I1 = IRANK(I)
         IRANK(I) = K
         K = K + 1
         I = I1
         IF (I.GT.0) GO TO 300
      END IF
C
      IF (IERR.NE.0) THEN
         IFAIL = P01ABF(IFAIL,IERR,SRNAME,2,P01REC)
      ELSE
         IFAIL = 0
      END IF
      RETURN
C
99999 FORMAT (' ** On entry, one or more of the following parameter va',
     *  'lues is illegal',/'    M1 =',I16,'  M2 =',I16)
      END
