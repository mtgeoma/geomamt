      SUBROUTINE M01ZBF(IPERM,M1,M2,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     M01ZBF CHECKS THE VALIDITY OF A PERMUTATION.
C
C     M01ZBF CHECKS THAT ELEMENTS M1 TO M2 OF IPERM CONTAIN A VALID
C     PERMUTATION OF THE INTEGERS M1 TO M2. THE CONTENTS OF IPERM ARE
C     UNCHANGED ON EXIT.
C
C     WRITTEN BY N.N.MACLAREN, UNIVERSITY OF CAMBRIDGE.
C     REVISED BY NAG CENTRAL OFFICE.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='M01ZBF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M1, M2
C     .. Array Arguments ..
      INTEGER           IPERM(M2)
C     .. Local Scalars ..
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
C       CHECK THE PARAMETERS.
C
      IF (M2.LT.1 .OR. M1.LT.1 .OR. M1.GT.M2) THEN
         IERR = 1
         WRITE (P01REC,FMT=99999) M1, M2
      ELSE
         IERR = 0
C
C        CHECK THAT IPERM IS WITHIN RANGE
C
         DO 20 I = M1, M2
            J = IPERM(I)
            IF ((J.LT.M1) .OR. (J.GT.M2)) GO TO 100
            IF (I.NE.J) IPERM(I) = -J
   20    CONTINUE
C
C        CHECK THAT NO VALUE IS REPEATED
C
         DO 60 I = M1, M2
            K = -IPERM(I)
            IF (K.GE.0) THEN
               J = I
   40          IPERM(J) = K
               J = K
               K = -IPERM(J)
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
C     RESTORE IPERM
C
      DO 160 I = M1, M2
         IPERM(I) = ABS(IPERM(I))
  160 CONTINUE
      GO TO 80
C
99999 FORMAT (' ** On entry, one or more of the following parameter va',
     *  'lues is illegal',/'    M1 =',I16,'  M2 =',I16)
99998 FORMAT (' ** IPERM(M1:M2) does not contain a permutation of the ',
     *  'integers M1 to M2')
99997 FORMAT ('    IPERM(',I6,') contains an out-of-range value =',I16)
99996 FORMAT ('    IPERM contains a duplicate value =',I16)
      END
