      INTEGER FUNCTION G02EAZ(M,N,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     BASED ON TOMS ALGORTIHM 515 WHICH WAS BASED ON
C     ACM ALGORITHM 160
C
C     CALCULATES THE NUMBER OF COMBINATIONS OF
C     M OUT OF N
C
C     .. Parameters ..
      CHARACTER*6             SRNAME
      PARAMETER               (SRNAME='G02EAZ')
C     .. Scalar Arguments ..
      INTEGER                 IFAIL, M, N
C     .. Local Scalars ..
      DOUBLE PRECISION        X
      INTEGER                 I, IERROR, MAX, NREC, P, Q, R
C     .. Local Arrays ..
      CHARACTER*80            P01REC(2)
C     .. External Functions ..
      INTEGER                 P01ABF, X02BBF
      EXTERNAL                P01ABF, X02BBF
C     .. Executable Statements ..
      G02EAZ = 0
      IF (M.GT.0 .AND. N.GE.M) THEN
         NREC = 1
         IERROR = 0
         MAX = X02BBF(X)
         Q = M
         P = N - Q
         IF (Q.LT.P) THEN
            P = Q
            Q = N - P
         END IF
         R = Q + 1
         IF (P.EQ.0) R = 1
         IF (P.GE.2) THEN
            DO 20 I = 2, P
               IF (R.GT.MAX/(Q+I)) THEN
                  GO TO 40
C
               ELSE
                  R = ((Q+I)*R)/I
               END IF
   20       CONTINUE
            GO TO 60
C
   40       IERROR = 2
            WRITE (P01REC,FMT=99998)
            R = MAX
         END IF
   60    G02EAZ = R
      ELSE
         IERROR = 1
         NREC = 2
         WRITE (P01REC,FMT=99999) M, N
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (' ** On entry values of M and N are invalid',/'  M = ',
     *       I16,' N = ',I16)
99998 FORMAT (' ** Overflow would occur in calculation')
      END
