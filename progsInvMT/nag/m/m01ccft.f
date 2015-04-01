      SUBROUTINE M01CCF(CH,M1,M2,L1,L2,ORDER,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     M01CCF RE-ARRANGES A VECTOR OF CHARACTER DATA SO THAT A
C     SPECIFIED SUBSTRING IS IN ASCII OR REVERSE ASCII ORDER.
C
C     M01CCF IS BASED ON SINGLETON'S IMPLEMENTATION OF THE
C     'MEDIAN-OF-THREE' QUICKSORT ALGORITHM, BUT WITH TWO
C     ADDITIONAL MODIFICATIONS. FIRST, SMALL SUBFILES ARE
C     SORTED BY AN INSERTION SORT ON A SEPARATE FINAL PASS.
C     SECOND, IF A LARGER OF THEM IS FLAGGED FOR SPECIAL
C     TREATMENT: BEFORE IT IS PARTITIONED, ITS END-POINTS
C     ARE SWAPPED WITH TWO RANDOM POINTS WITHIN IT; THIS
C     MAKES THE WORST CASE BEHAVIOUR EXTREMELY UNLIKELY.
C
C     ONLY THE SUBSTRING (L1:L2) OF EACH ELEMENT OF THE ARRAY
C     CH IS USED TO DETERMINE THE SORTED ORDER, BUT THE ENTIRE
C     ELEMENTS ARE RE-ARRANGED INTO SORTED ORDER.
C
C     THE MAXIMUM LENGTH OF A SMALL SUBFILE IS DEFINED BY THE
C     VARIABLE MINQIK, SET TO 15.
C
C     THE ROUTINE ASSUMES THAT THE NUMBER OF ELEMENTS TO BE
C     SORTED DOES NOT EXCEED MINQIK*2**MAXSTK.
C
C     THE MAXIMUM PERMITTED LENGTH OF EACH ELEMENT OF THE ARRAY CH
C     IS DEFINED BY THE PARAMETER MAXLCH. THIS RESTRICTION IS
C     IMPOSED BY THE NEED TO SPECIFY A LENGTH FOR THE INTERNAL
C     CHARACTER VARIABLE A.
C
C     WRITTEN BY NAG CENTRAL OFFICE.
C
C     .. Parameters ..
      INTEGER           MAXSTK
      PARAMETER         (MAXSTK=40)
      INTEGER           MINQIK
      PARAMETER         (MINQIK=15)
      INTEGER           MAXLCH
      PARAMETER         (MAXLCH=255)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='M01CCF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, L1, L2, M1, M2
      CHARACTER*1       ORDER
C     .. Array Arguments ..
      CHARACTER*(*)     CH(M2)
C     .. Local Scalars ..
      INTEGER           I, I1, I2, IERR, IR1, IR2, IR3, ISTK, J, J1, J2,
     *                  K, LCH, LENG, LL1, LL2, NREC, RAND
      CHARACTER*(MAXLCH) A
C     .. Local Arrays ..
      INTEGER           IHIGH(MAXSTK), ILOW(MAXSTK)
      CHARACTER*80      P01REC(4)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         MOD, LEN, LGT, LLT
C     .. Save statement ..
      SAVE              IR1, IR2, IR3
C     .. Data statements ..
      DATA              IR1, IR2, IR3/15223, 17795, 28707/
C     .. Executable Statements ..
C
C       CHECK THE PARAMETERS AND DECIDE IF QUICKSORT IS NEEDED.
C
      LCH = LEN(CH(1))
      IF (M2.LT.1 .OR. M1.LT.1 .OR. M1.GT.M2 .OR. L2.LT.1 .OR. L1.LT.
     *    1 .OR. L1.GT.L2 .OR. L2.GT.LCH) THEN
         IERR = 1
         WRITE (P01REC,FMT=99999) M1, M2, L1, L2
         NREC = 3
         IF (L2.GT.LCH) THEN
            WRITE (P01REC(4),FMT=99996) LCH
            NREC = 4
         END IF
      ELSE IF (ORDER.NE.'A' .AND. ORDER.NE.'a' .AND. ORDER.NE.'R' .AND.
     *         ORDER.NE.'r') THEN
         IERR = 2
         WRITE (P01REC,FMT=99998) ORDER
         NREC = 1
      ELSE IF (LCH.GT.MAXLCH) THEN
         IERR = 3
         WRITE (P01REC,FMT=99997) MAXLCH, LCH
         NREC = 1
      ELSE IF (M1.EQ.M2) THEN
         IERR = 0
      ELSE
         IERR = 0
         LENG = M2 - M1 + 1
         LL1 = L1
         LL2 = L2
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
            RAND = MOD(IR1/30269+IR2/30307+IR3/30323,1)
            K = I + RAND*(J-I)
            A(1:LCH) = CH(I)
            CH(I) = CH(K)
            CH(K) = A(1:LCH)
            K = I + J - K
            A(1:LCH) = CH(K)
            CH(K) = CH(J)
            CH(J) = A(1:LCH)
         END IF
C
C           CALCULATE A MEDIAN BY SINGLETONS METHOD.
C
         K = (I+J)/2
         IF (LGT(CH(I)(LL1:LL2),CH(J)(LL1:LL2))) THEN
            A(1:LCH) = CH(I)
            CH(I) = CH(J)
            CH(J) = A(1:LCH)
         END IF
         A(1:LCH) = CH(K)
         IF (LLT(A(LL1:LL2),CH(I)(LL1:LL2))) THEN
            CH(K) = CH(I)
            CH(I) = A(1:LCH)
            A(1:LCH) = CH(K)
         ELSE IF (LGT(A(LL1:LL2),CH(J)(LL1:LL2))) THEN
            CH(K) = CH(J)
            CH(J) = A(1:LCH)
            A(1:LCH) = CH(K)
         END IF
C
C           SPLIT THE VECTOR INTO TWO ASCENDING PARTS.  THIS IS WHERE
C           THE TIME IS SPENT.
C
         I1 = I
         J1 = J
   40    I1 = I1 + 1
         IF (LLT(CH(I1)(LL1:LL2),A(LL1:LL2))) GO TO 40
   60    J1 = J1 - 1
         IF (LGT(CH(J1)(LL1:LL2),A(LL1:LL2))) GO TO 60
         IF (I1.GE.J1) GO TO 80
         A(1:LCH) = CH(I1)
         CH(I1) = CH(J1)
         CH(J1) = A(1:LCH)
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
            A(1:LCH) = CH(I)
            J = I - 1
            IF (LLT(A(LL1:LL2),CH(J)(LL1:LL2))) THEN
  120          CH(J+1) = CH(J)
               J = J - 1
               IF (J.GE.M1) THEN
                  IF (LLT(A(LL1:LL2),CH(J)(LL1:LL2))) GO TO 120
               END IF
               CH(J+1) = A(1:LCH)
            END IF
  140    CONTINUE
C
C           REVERSE THE ORDER IF NECESSARY AND RETURN.
C
         IF (ORDER.EQ.'R' .OR. ORDER.EQ.'r') THEN
            DO 160 I = M1, (M1+M2-1)/2
               I1 = M1 + M2 - I
               A(1:LCH) = CH(I)
               CH(I) = CH(I1)
               CH(I1) = A(1:LCH)
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
     *  'lues is illegal',/'    M1 =',I16,'  M2 =',I16,/'    L1 =',I16,
     *  '  L2 =',I16)
99998 FORMAT (' ** On entry, ORDER has an illegal value: ORDER = ',A1)
99997 FORMAT (' ** LEN(CH(1)) is greater than ',I3,': LEN(CH(1)) =',I16)
99996 FORMAT ('    L2 must not be greater than LEN(CH(1)): LEN(CH(1)) ='
     *  ,I6)
      END
