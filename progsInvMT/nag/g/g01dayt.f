      DOUBLE PRECISION FUNCTION G01DAY(IS,N,IERROR)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Computes Normal Scores
C
C     A single call single iteration version of G01DAF
C
C     .. Scalar Arguments ..
      INTEGER                          IERROR, IS, N
C     .. Local Scalars ..
      DOUBLE PRECISION                 CRLN, ERR, RINC, SCORE, TWOPI
      INTEGER                          I, J, M
      LOGICAL                          SIGN
C     .. External Functions ..
      DOUBLE PRECISION                 X01AAF
      EXTERNAL                         X01AAF
C     .. External Subroutines ..
      EXTERNAL                         G01DAX
C     .. Intrinsic Functions ..
      INTRINSIC                        LOG, DBLE
C     .. Executable Statements ..
      IERROR = 0
      M = N/2
      SIGN = .FALSE.
C
C      switch for IS gt M/2
C
      IF (IS.LE.M) THEN
         I = IS
      ELSE IF (IS.EQ.M+1 .AND. 2*M.LT.N) THEN
         G01DAY = 0.0D0
         GO TO 40
      ELSE
         SIGN = .TRUE.
         I = N - IS + 1
      END IF
C
C     set up constants
C
      TWOPI = 2.0D0*X01AAF(TWOPI)
      CRLN = LOG(DBLE(N)) - 0.5D0*LOG(TWOPI)
      DO 20 J = 1, I - 1
         CRLN = CRLN + LOG(DBLE(N-J)) - LOG(DBLE(J))
   20 CONTINUE
      IF (N.GE.300) THEN
         RINC = 0.01D0
      ELSE IF (N.GE.1100) THEN
         RINC = 0.005D0
      ELSE
         RINC = 0.02D0
      END IF
C
C     start calculations (single iteration)
C
      CALL G01DAX(I,N,CRLN,RINC,SCORE,ERR)
      IF (SIGN) THEN
         G01DAY = -SCORE
      ELSE
         G01DAY = SCORE
      END IF
   40 RETURN
      END
