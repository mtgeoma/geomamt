      SUBROUTINE G11SBF(N2,N,S,X,NRX,RL,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G11SBF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, N, N2, NRX, S
C     .. Array Arguments ..
      INTEGER           RL(N)
      LOGICAL           X(NRX,N2)
C     .. Local Scalars ..
      INTEGER           I, IERROR, J, K
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
C
C     TEST FOR ERRORS IN INPUT DATA
C
      IF (N2.LT.3) GO TO 20
      IF (N.LT.7) GO TO 20
      IF (NRX.LT.N) GO TO 20
C
C     IF WE HAVE GOT TO THIS POINT THEN THE INPUT DATA
C     IS FREE FROM ERROR AND WE DO NOT NEED TO RESET IFAIL
C
      GO TO 40
C
C
   20 IERROR = 1
      WRITE (P01REC,FMT=99999) N2, N, NRX
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,2,P01REC)
      RETURN
C
C     COUNT THE OBSERVED FREQUENCIES OF THE S DIFFERENT SCORE PATTERNS
C     ,ASSIGNING TO RL,AND THEN RE-ORDER THE ROWS OF X SUCH THAT
C     THE FIRST S ROWS OF X DISPLAY THE S DIFFERENT SCORE PATTERNS
C
   40 S = 0
C
C     INITIALISE SCORE COUNTS IN RL TO ONE
C
      DO 60 I = 1, N
         RL(I) = 1
   60 CONTINUE
C
      DO 120 I = 1, N - 1
C
C        CHECK WHETHER I TH ROW HAS ALREADY BEEN COUNTED
C        (RL(I) SET EQUAL TO ZERO) AND IF SO LOOK AT I+1 TH
C
         IF (RL(I).EQ.0) GO TO 120
C
C        I TH ROW OF X IS A NEW SCORE PATTERN
C
         S = S + 1
C
C        COUNT THE NUMBER OF TIMES THE I TH SCORE PATTERN IS
C        REPLICATED IN THE ROWS BEYOND THE I TH
C
         DO 100 K = I + 1, N
C
C           IF THE K TH ROW HAS ALREADY BEEN COUNTED PREVIOUSLY
C           WE PROCEED TO THE K+1 TH ROW
C
            IF (RL(K).EQ.0) GO TO 100
C
            DO 80 J = 1, N2
               IF (X(K,J) .NEQV. X(I,J)) GO TO 100
   80       CONTINUE
C
C           K TH SCORE PATTERN IS IDENTICAL TO THE I TH SO INCREMENT
C           COUNT ON I TH SCORE PATTERN AND REGISTER THAT THE K TH ROW
C           HAS NOW BEEN COUNTED I.E. SET THE K TH COMPONENT OF
C           RL TO ZERO
C
            RL(I) = RL(I) + 1
            RL(K) = 0
C
  100    CONTINUE
C
  120 CONTINUE
C
      IF (RL(N).EQ.1) S = S + 1
C
C     SHIFT THE S DIFFERENT SCORE PATTERNS UP TO THE FIRST S ROWS OF X
C     AND SORT THE ELEMENTS OF RL ACCORDINGLY. IF ALL SCORE
C     PATTERNS ARE DIFFERENT NO SORTING NEEDS TO BE DONE.
C
      IF (S.EQ.N) GO TO 180
      K = 2
      DO 160 I = 2, N
         IF (RL(I).EQ.0) GO TO 160
         RL(K) = RL(I)
         DO 140 J = 1, N2
            X(K,J) = X(I,J)
  140    CONTINUE
         K = K + 1
  160 CONTINUE
C
  180 IFAIL = 0
      RETURN
C
99999 FORMAT (' ** ON ENTRY, ONE OR MORE OF THE FOLLOWING PARAMETER VA',
     *  'LUES IS ILLEGAL',/'   IP = ',I16,'   N = ',I16,'   NRX = ',I16)
      END
