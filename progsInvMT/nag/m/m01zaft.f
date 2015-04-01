      SUBROUTINE M01ZAF(IPERM,M1,M2,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     M01ZAF INVERTS A PERMUTATION, AND HENCE CONVERTS A RANK VECTOR TO
C     AN INDEX VECTOR, OR VICE VERSA.
C
C     WRITTEN BY N.M.MACLAREN, UNIVERSITY OF CAMBRIDGE.
C     REVISED BY NAG CENTRAL OFFICE.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='M01ZAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M1, M2
C     .. Array Arguments ..
      INTEGER           IPERM(M2)
C     .. Local Scalars ..
      INTEGER           I, IERR, J, K, L
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
C
C       CHECK THE PARAMETERS AND MODIFY IPERM.
C
      IF (M2.LT.1 .OR. M1.LT.1 .OR. M1.GT.M2) THEN
         IERR = 1
         WRITE (P01REC,FMT=99999) M1, M2
      ELSE
         IERR = 0
         DO 20 I = M1, M2
            J = IPERM(I)
            IF ((J.LT.M1) .OR. (J.GT.M2)) GO TO 100
            IF (I.NE.J) IPERM(I) = -J
   20    CONTINUE
C
C        REVERSE EACH NON-TRIVIAL CYCLE.
C
         DO 60 I = M1, M2
            K = -IPERM(I)
            IF (K.GE.0) THEN
               J = I
   40          L = -IPERM(K)
               IF (L.GT.0) THEN
                  IPERM(K) = J
                  J = K
                  K = L
                  GO TO 40
               END IF
               IF (IPERM(I).LT.0) GO TO 120
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
      WRITE (P01REC(2),FMT=99996) K
  140 WRITE (P01REC(1),FMT=99998)
      GO TO 80
C
99999 FORMAT (' ** On entry, one or more of the following parameter va',
     *  'lues is illegal',/'    M1 =',I16,'  M2 =',I16)
99998 FORMAT (' ** IPERM(M1:M2) does not contain a permutation of the ',
     *  'integers M1 to M2')
99997 FORMAT ('    IPERM(',I6,') contains an out-of-range value',I16)
99996 FORMAT ('    IPERM contains a repeated value',I16)
      END
