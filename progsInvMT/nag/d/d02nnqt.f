      SUBROUTINE D02NNQ(MSG,IERT,NI,I1,I2,NR,R1,R2)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-713 (DEC 1989).
C
C     OLD NAME SERROR
C
C-----------------------------------------------------------------------
C  ERROR HANDLING ROUTINE FOR THE SPRINT INTEGRATION PACKAGE. THIS
C  ROUTINE IS A FORTRAN77 IMPROVED VERSION OF THE ROUTINE USED IN LSODI
C  AND MAKES USE OF CHARACTER HANDLING FACILITIES.
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  R1, R2
      INTEGER           I1, I2, IERT, NI, NR
      CHARACTER*(*)     MSG
C     .. Scalars in Common ..
      INTEGER           ITRACE, NERR
C     .. Local Scalars ..
      INTEGER           I, IL, IT, J, K, KP1, LWORD, NERR1, NERR2
      CHARACTER*80      REC
      CHARACTER*(240)   MSG1
C     .. Local Arrays ..
      CHARACTER*(60)    MSGOUT(5)
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN, LEN
C     .. Common blocks ..
      COMMON            /AD02NM/ITRACE, NERR
C     .. Save statement ..
      SAVE              /AD02NM/
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C ALL ARGUMENTS ARE INPUT ARGUMENTS.
C
C MSG    = THE MESSAGE IN CHARACTER FORMAT
C NI     = NUMBER OF INTEGERS (0, 1, OR 2) TO BE PRINTED WITH MESSAGE.
C I1,I2  = INTEGERS TO BE PRINTED, DEPENDING ON NI.
C NR     = NUMBER OF REALS (0, 1, OR 2) TO BE PRINTED WITH MESSAGE.
C R1,R2  = REALS TO BE PRINTED, DEPENDING ON NR.
C-----------------------------------------------------------------------
C
      IF (ITRACE.LT.0 .OR. ITRACE.EQ.0 .AND. IERT.EQ.2) RETURN
      CALL X04AAF(0,NERR1)
      CALL X04ABF(0,NERR2)
C
      IL = LEN(MSG)
C
C    SET MSG1 BLANK AND GET RID OF UNNECESSARY SPACES IN ERROR MESSAGE
C
      J = 1
      IT = MIN(IL,240)
      DO 20 I = 1, 10
         MSG1(J:) = '                        '
         J = J + 24
   20 CONTINUE
      K = 0
      J = 0
      DO 40 I = 1, IT
         IF (MSG(I:I).EQ.' ') THEN
            K = K + 1
            IF (K.GT.1) GO TO 40
         ELSE
            K = 0
         END IF
         J = J + 1
         MSG1(J:J) = MSG(I:I)
   40 CONTINUE
      IL = J
C
C     FORMAT THE MESSAGE NOW STORED IN MSG1
C
      I = 1
      J = 0
   60 LWORD = I + 60
      J = J + 1
      IF (IL.LE.LWORD) THEN
         MSGOUT(J) = MSG1(I:IL)
         GO TO 100
      ELSE IF (MSG1(LWORD:LWORD).EQ.' ') THEN
         MSGOUT(J) = MSG1(I:LWORD)
      ELSE
         DO 80 K = LWORD, I, -1
            IF (MSG1(K:K).EQ.' ') GO TO 90
   80    CONTINUE
         K = LWORD
   90    MSGOUT(J) = MSG1(I:K)
         LWORD = K
      END IF
      I = LWORD
      GO TO 60
C
C  OUTPUT THE ERROR MESSAGE
C
  100 DO 120 I = 1, J
         IF (IERT.EQ.1) THEN
            CALL X04BAF(NERR1,MSGOUT(I))
         ELSE
            CALL X04BAF(NERR2,MSGOUT(I))
         END IF
  120 CONTINUE
C
C  PRINT THE INTEGERS AND REALS IN THE ERROR MESSAGE (IF ANY)
C
      IF (NI.EQ.1) THEN
         WRITE (REC,FMT=99999) I1
         IF (IERT.EQ.1) THEN
            CALL X04BAF(NERR1,REC)
         ELSE
            CALL X04BAF(NERR2,REC)
         END IF
      ELSE IF (NI.EQ.2) THEN
         WRITE (REC,FMT=99998) I1, I2
         IF (IERT.EQ.1) THEN
            CALL X04BAF(NERR1,REC)
         ELSE
            CALL X04BAF(NERR2,REC)
         END IF
      END IF
      IF (NR.EQ.1) THEN
         WRITE (REC,FMT=99997) R1
         IF (IERT.EQ.1) THEN
            CALL X04BAF(NERR1,REC)
         ELSE
            CALL X04BAF(NERR2,REC)
         END IF
      ELSE IF (NR.EQ.2) THEN
         WRITE (REC,FMT=99996) R1, R2
         IF (IERT.EQ.1) THEN
            CALL X04BAF(NERR1,REC)
         ELSE
            CALL X04BAF(NERR2,REC)
         END IF
      END IF
      RETURN
C
99999 FORMAT (' IN ABOVE MESSAGE I1 =',I10)
99998 FORMAT (' IN ABOVE MESSAGE I1 =',I10,'   I2 =',I10)
99997 FORMAT (' IN ABOVE MESSAGE R1 =',D15.7)
99996 FORMAT (' IN ABOVE MESSAGE R1 =',D15.7,'   R2 =',D15.7)
      END
