      SUBROUTINE M01EAF(RV,M1,M2,IRANK,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     M01EAF RE-ARRANGES A VECTOR OF REAL NUMBERS INTO
C     THE ORDER SPECIFIED BY A VECTOR OF RANKS.
C
C     M01EAF IS DESIGNED TO BE USED TYPICALLY IN CONJUNCTION
C     WITH THE M01D- RANKING ROUTINES. AFTER ONE OF THE M01D-
C     ROUTINES HAS BEEN CALLED TO DETERMINE A VECTOR OF RANKS,
C     M01EAF CAN BE CALLED TO RE-ARRANGE A VECTOR OF REAL
C     NUMBERS INTO THE RANK ORDER. IF THE VECTOR OF RANKS HAS
C     BEEN GENERATED IN SOME OTHER WAY, THEN M01ZBF SHOULD BE
C     CALLED TO CHECK ITS VALIDITY BEFORE M01EAF IS CALLED.
C
C     WRITTEN BY N.M.MACLAREN, UNIVERSITY OF CAMBRIDGE.
C     REVISED BY NAG CENTRAL OFFICE.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='M01EAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M1, M2
C     .. Array Arguments ..
      DOUBLE PRECISION  RV(M2)
      INTEGER           IRANK(M2)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B
      INTEGER           I, IERR, J, K
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
C
C       CHECK THE PARAMETERS AND MODIFY IRANK.
C
      IF (M2.LT.1 .OR. M1.LT.1 .OR. M1.GT.M2) THEN
         IERR = 1
         WRITE (P01REC,FMT=99999) M1, M2
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
               A = RV(I)
   40          IRANK(J) = K
               B = RV(K)
               RV(K) = A
               J = K
               A = B
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
         IFAIL = P01ABF(IFAIL,IERR,SRNAME,2,P01REC)
      ELSE
         IFAIL = 0
      END IF
      RETURN
  100 IERR = 2
      WRITE (P01REC(2),FMT=99997) I, J
      GO TO 140
  120 IERR = 3
      WRITE (P01REC(2),FMT=99996) J
  140 WRITE (P01REC(1),FMT=99998)
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
99997 FORMAT ('    IRANK(',I6,') contains an out-of-range value',I16)
99996 FORMAT ('    IRANK contains a repeated value',I16)
      END
