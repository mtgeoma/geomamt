      SUBROUTINE G08EAY(UPPER,N,A,LDA,X,STAT,WK,IERROR)
C     MARK 16 REVISED. IER-1040 (JUN 1993).
C
C     COMPUTES S=X'.INV(A).X
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  STAT
      INTEGER           IERROR, LDA, N
      CHARACTER*1       UPPER
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,N), WK((N*N+3*N)/2), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS
      INTEGER           I, IFAULT, IJ, IL, J, N1
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          F01BQZ, DCOPY, DTPSV
C     .. Executable Statements ..
      IERROR = 0
      IF (N.LE.0 .OR. LDA.LT.N) THEN
         IERROR = 1
      ELSE IF (N.EQ.1) THEN
         IF (A(1,1).GT.0.0D0) THEN
            STAT = X(1)*X(1)/A(1,1)
         ELSE
            IERROR = 2
         END IF
      ELSE
         IF (UPPER.EQ.'U' .OR. UPPER.EQ.'u') THEN
            IJ = N
            DO 40 I = 1, N
               WK(I) = A(I,I)
               DO 20 J = 1, I - 1
                  IJ = IJ + 1
                  WK(IJ) = A(J,I)
   20          CONTINUE
   40       CONTINUE
         ELSE IF (UPPER.EQ.'L' .OR. UPPER.EQ.'l') THEN
            IJ = N
            DO 80 I = 1, N
               WK(I) = A(I,I)
               DO 60 J = 1, I - 1
                  IJ = IJ + 1
                  WK(IJ) = A(I,J)
   60          CONTINUE
   80       CONTINUE
         ELSE
            IERROR = 1
            GO TO 160
         END IF
C
         EPS = X02AJF()
         N1 = (N*N-N)/2
         IFAULT = 1
         CALL F01BQZ(N,EPS,WK(N+1),N1,WK,IFAULT)
         IF (IFAULT.NE.0) THEN
            IERROR = 2
         END IF
         IJ = N1 + 2*N
         IL = N + N1
         DO 120 I = N, 2, -1
            WK(IJ) = WK(I)
            IJ = IJ - 1
            DO 100 J = I - 1, 1, -1
               WK(IJ) = WK(IL)
               IL = IL - 1
               IJ = IJ - 1
  100       CONTINUE
  120    CONTINUE
         WK(N+1) = WK(1)
         CALL DCOPY(N,X,1,WK,1)
         CALL DTPSV('U','T','U',N,WK(N+1),WK,1)
         IJ = N
         STAT = 0.0D0
         DO 140 I = 1, N
            IJ = IJ + I
            STAT = (WK(I)**2)/WK(IJ) + STAT
  140    CONTINUE
      END IF
  160 RETURN
C
      END
