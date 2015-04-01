      SUBROUTINE G07BEY(RV,M1,M2,ORDER,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     G07BEY SORTS A VECTOR OF REAL NUMBERS INTO ASCENDING
C     OR DESCENDING ORDER USING ABSOLUTE VALUES. G07BEY IS
C     BASED ON M01CAF.
C
C     M01CAF IS BASED ON SINGLETON'S IMPLEMENTATION OF THE
C     'MEDIAN-OF-THREE' QUICKSORT ALGORITHM, BUT WITH TWO
C     ADDITIONAL MODIFICATIONS. FIRST, SMALL SUBFILES ARE
C     SORTED BY AN INSERTION SORT ON A SEPARATE FINAL PASS.
C     SECOND, IF A SUBFILE IS PARTITIONED INTO TWO VERY
C     UNBALANCED SUBFILES, THE LARGER OF THEM IS FLAGGED FOR
C     SPECIAL TREATMENT: BEFORE IT IS PARTITIONED, ITS END-
C     POINTS ARE SWAPPED WITH TWO RANDOM POINTS WITHIN IT;
C     THIS MAKES THE WORST CASE BEHAVIOUR EXTREMELY UNLIKELY.
C
C     THE MAXIMUM LENGTH OF A SMALL SUBFILE IS DEFINED BY THE
C     VARIABLE MINQIK, SET TO 15.
C
C     THE ROUTINE ASSUMES THAT THE NUMBER OF ELEMENTS TO BE
C     SORTED DOES NOT EXCEED MINQIK*2**MAXSTK.
C
C     .. Parameters ..
      INTEGER           MAXSTK
      PARAMETER         (MAXSTK=40)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G07BEY')
      INTEGER           MINQIK
      PARAMETER         (MINQIK=15)
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M1, M2
      CHARACTER         ORDER
C     .. Array Arguments ..
      DOUBLE PRECISION  RV(M2)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, RAND, X
      INTEGER           I, I1, I2, IERR, IR1, IR2, IR3, ISTK, J, J1, J2,
     *                  K, LENG, NREC
C     .. Local Arrays ..
      INTEGER           IHIGH(MAXSTK), ILOW(MAXSTK)
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MOD, DBLE
C     .. Save statement ..
      SAVE              IR1, IR2, IR3
C     .. Data statements ..
      DATA              IR1, IR2, IR3/15223, 17795, 28707/
C     .. Executable Statements ..
C
C       CHECK THE PARAMETERS AND DECIDE IF QUICKSORT IS NEEDED.
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
         IERR = 0
      ELSE
         IERR = 0
         LENG = M2 - M1 + 1
         IF (LENG.LE.MINQIK) GO TO 100
C
C           INITIALISE AND START QUICKSORT ON THE WHOLE VECTOR.
C
         ISTK = 0
         I = M1
         J = M2
C
C           IF THE PREVIOUS PASS WAS BAD, CHANGE THE END VALUES AT
C           RANDOM.
C
   20    IF (I.LT.0) THEN
            I = -I
            IR1 = 171*MOD(IR1,177) - 2*(IR1/177)
            IR2 = 172*MOD(IR2,176) - 35*(IR2/176)
            IR3 = 170*MOD(IR3,178) - 63*(IR3/178)
            IF (IR1.LT.0) IR1 = IR1 + 30269
            IF (IR2.LT.0) IR2 = IR2 + 30307
            IF (IR3.LT.0) IR3 = IR3 + 30323
            RAND = MOD(DBLE(IR1)/30269.0D0+DBLE(IR2)/30307.0D0+DBLE(IR3)
     *             /30323.0D0,1.0D0)
            K = I + RAND*(J-I)
            X = RV(I)
            RV(I) = RV(K)
            RV(K) = X
            K = I + J - K
            X = RV(K)
            RV(K) = RV(J)
            RV(J) = X
         END IF
C
C           CALCULATE A MEDIAN BY SINGLETONS METHOD.
C
         K = (I+J)/2
         IF (ABS(RV(I)).GT.ABS(RV(J))) THEN
            X = RV(I)
            RV(I) = RV(J)
            RV(J) = X
         END IF
         A = RV(K)
         IF (ABS(A).LT.ABS(RV(I))) THEN
            RV(K) = RV(I)
            RV(I) = A
            A = RV(K)
         ELSE IF (ABS(A).GT.ABS(RV(J))) THEN
            RV(K) = RV(J)
            RV(J) = A
            A = RV(K)
         END IF
C
C           SPLIT THE VECTOR INTO TWO ASCENDING PARTS.  THIS IS WHERE
C           THE TIME IS SPENT.
C
         I1 = I
         J1 = J
   40    I1 = I1 + 1
         IF (ABS(RV(I1)).LT.ABS(A)) GO TO 40
   60    J1 = J1 - 1
         IF (ABS(RV(J1)).GT.ABS(A)) GO TO 60
         IF (I1.GE.J1) GO TO 80
         X = RV(I1)
         RV(I1) = RV(J1)
         RV(J1) = X
         GO TO 40
C
C           STACK ONE SUBFILE, IF APPROPRIATE, AND CARRY ON.
C
   80    I2 = I1 - I
         J2 = J - J1
         IF (J2.LE.I2) THEN
            IF (I2.LE.MINQIK) THEN
               IF (ISTK.LE.0) GO TO 100
               I = ILOW(ISTK)
               J = IHIGH(ISTK)
               ISTK = ISTK - 1
            ELSE
C
C                 TEST FOR VERY UNBALANCED SUBFILES
C                 ( THE DETAILS OF THE TEST ARE FAIRLY ARBITRARY.)
C
               IF (5*(J2+5).LT.I2) I = -I
               IF (J2.LE.MINQIK) THEN
                  J = I1 - 1
               ELSE
                  ISTK = ISTK + 1
                  ILOW(ISTK) = I
                  IHIGH(ISTK) = I1 - 1
                  I = J1 + 1
               END IF
            END IF
         ELSE
C
C              DEAL WITH THE CASE WHEN THE SECOND PART IS LARGER.
C
            IF (J2.LE.MINQIK) THEN
               IF (ISTK.LE.0) GO TO 100
               I = ILOW(ISTK)
               J = IHIGH(ISTK)
               ISTK = ISTK - 1
            ELSE
C
C                 TEST FOR VERY UNBALANCED SUBFILES
C                 ( THE DETAILS OF THE TEST ARE FAIRLY ARBITRARY.)
C
               IF (5*(I2+5).LT.J2) J1 = -(J1+2)
               IF (I2.LE.MINQIK) THEN
                  I = J1 + 1
               ELSE
                  ISTK = ISTK + 1
                  ILOW(ISTK) = J1 + 1
                  IHIGH(ISTK) = J
                  J = I1 - 1
               END IF
            END IF
         END IF
         GO TO 20
C
C           TIDY UP AND DO AN ASCENDING INSERTION SORT.
C
  100    DO 140 I = M1 + 1, M2
            A = RV(I)
            J = I - 1
            IF (ABS(A).LT.ABS(RV(J))) THEN
  120          RV(J+1) = RV(J)
               J = J - 1
               IF (J.GE.M1) THEN
                  IF (ABS(A).LT.ABS(RV(J))) GO TO 120
               END IF
               RV(J+1) = A
            END IF
  140    CONTINUE
C
C           REVERSE THE ORDER IF NECESSARY AND RETURN.
C
         IF ((ORDER.EQ.'D') .OR. (ORDER.EQ.'d')) THEN
            DO 160 I = M1, (M1+M2-1)/2
               I1 = M1 + M2 - I
               X = RV(I)
               RV(I) = RV(I1)
               RV(I1) = X
  160       CONTINUE
         END IF
C
      END IF
      IF (IERR.NE.0) THEN
         IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,P01REC)
      ELSE
         IFAIL = 0
      END IF
      RETURN
C
99999 FORMAT (' ** On entry, one or more of the following parameter va',
     *       'lues is illegal',/'    M1 =',I16,'  M2 =',I16)
99998 FORMAT (' ** On entry, ORDER has an illegal value: ORDER = ',A1)
      END
