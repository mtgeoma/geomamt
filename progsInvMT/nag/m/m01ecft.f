      SUBROUTINE M01ECF(CH,M1,M2,IRANK,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     M01ECF RE-ARRANGES A VECTOR OF CHARACTER DATA INTO
C     THE ORDER SPECIFIED BY A VECTOR OF RANKS.
C
C     M01ECF IS DESIGNED TO BE USED TYPICALLY IN CONJUNCTION
C     WITH THE M01D- RANKING ROUTINES. AFTER ONE OF THE M01D-
C     ROUTINES HAS BEEN CALLED TO DETERMINE A VECTOR OF RANKS,
C     M01ECF CAN BE CALLED TO RE-ARRANGE A VECTOR OF CHARACTER
C     DATA INTO THE RANK ORDER. IF THE VECTOR OF RANKS HAS
C     BEEN GENERATED IN SOME OTHER WAY, THEN M01ZBF SHOULD BE
C     CALLED TO CHECK ITS VALIDITY BEFORE M01ECF IS CALLED.
C
C     THE MAXIMUM PERMITTED LENGTH OF EACH ELEMENT OF THE ARRAY CH
C     IS DEFINED BY THE PARAMETER MAXLCH. THIS RESTRICTION IS
C     IMPOSED BY THE NEED TO SPECIFY A LENGTH FOR THE INTERNAL
C     CHARACTER VARIABLES A AND B.
C
C     WRITTEN BY NAG CENTRAL OFFICE
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='M01ECF')
      INTEGER           MAXLCH
      PARAMETER         (MAXLCH=255)
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M1, M2
C     .. Array Arguments ..
      INTEGER           IRANK(M2)
      CHARACTER*(*)     CH(M2)
C     .. Local Scalars ..
      INTEGER           I, IERR, J, K, LCH, NREC
      CHARACTER*(MAXLCH) A, B
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, LEN
C     .. Executable Statements ..
C
C       CHECK THE PARAMETERS AND MODIFY IRANK.
C
      LCH = LEN(CH(1))
      IF (M2.LT.1 .OR. M1.LT.1 .OR. M1.GT.M2) THEN
         IERR = 1
         WRITE (P01REC,FMT=99999) M1, M2
         NREC = 2
      ELSE IF (LCH.GT.MAXLCH) THEN
         IERR = 2
         WRITE (P01REC,FMT=99995) MAXLCH, LCH
         NREC = 1
      ELSE
         IERR = 0
         DO 20 I = M1, M2
            J = IRANK(I)
            IF ((J.LT.M1) .OR. (J.GT.M2)) GO TO 100
            IF (I.NE.J) IRANK(I) = -J
   20    CONTINUE
C
C        MOVE EACH NON-TRIVIAL CYCLE ROUND.
C
         DO 60 I = M1, M2
            K = -IRANK(I)
            IF (K.GE.0) THEN
               J = I
               A(1:LCH) = CH(I)
   40          IRANK(J) = K
               B(1:LCH) = CH(K)
               CH(K) = A(1:LCH)
               J = K
               A(1:LCH) = B(1:LCH)
               K = -IRANK(J)
               IF (K.GT.0) GO TO 40
               IF (J.NE.I) GO TO 120
            END IF
   60    CONTINUE
      END IF
C
C       RETURN
C
   80 IF (IERR.NE.0) THEN
         IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,P01REC)
      ELSE
         IFAIL = 0
      END IF
      RETURN
  100 IERR = 3
      WRITE (P01REC(2),FMT=99997) I, J
      GO TO 140
  120 IERR = 4
      WRITE (P01REC(2),FMT=99996) J
  140 WRITE (P01REC(1),FMT=99998)
      NREC = 2
C
C     RESTORE IRANK
C
      DO 160 J = M1, M2
         IRANK(J) = ABS(IRANK(J))
  160 CONTINUE
      GO TO 80
C
99999 FORMAT (' ** On entry, one or more of the following parameter va',
     *  'lues is illegal',/'    M1 =',I16,'  M2 =',I16)
99998 FORMAT (' ** IRANK(M1:M2) does not contain a permutation of the ',
     *  'integers M1 to M2')
99997 FORMAT ('    IRANK(',I6,') contains an out=of-range value',I16)
99996 FORMAT ('    IRANK contains a repeated value',I16)
99995 FORMAT (' ** LEN(CH(1)) is greater than ',I3,': LEN(CH(1)) =',I16)
      END
