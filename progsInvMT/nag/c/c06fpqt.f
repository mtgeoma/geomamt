      SUBROUTINE C06FPQ(M,N,INIT,TRIG,Q,NQ,IERROR)
CVD$R NOVECTOR
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER           IERROR, M, N, NQ
      CHARACTER*1       INIT
C     .. Array Arguments ..
      DOUBLE PRECISION  TRIG(2*N)
      INTEGER           Q(30)
C     .. Local Scalars ..
      INTEGER           NCHECK
C     .. External Subroutines ..
      EXTERNAL          C06FPY, C06FPZ
C     .. Intrinsic Functions ..
      INTRINSIC         NINT
C     .. Save statement ..
      SAVE              NCHECK
C     .. Data statements ..
      DATA              NCHECK/-1/
C     .. Executable Statements ..
      IERROR = 0
C
      IF (M.LT.1) THEN
         IERROR = 1
         RETURN
      ELSE IF (N.LT.1) THEN
         IERROR = 2
         RETURN
      END IF
      IF (INIT.NE.'I' .AND. INIT.NE.'i' .AND. INIT.NE.'S' .AND. INIT.NE.
     *    's' .AND. INIT.NE.'R' .AND. INIT.NE.'r') THEN
         IERROR = 3
         RETURN
      END IF
      IF (INIT.EQ.'S' .OR. INIT.EQ.'s') THEN
         IF (NCHECK.EQ.-1) THEN
            IERROR = 4
            RETURN
         ELSE IF (NINT(TRIG(N)).NE.N .OR. NINT(TRIG(2*N)).NE.N) THEN
            IERROR = 5
            RETURN
         END IF
      END IF
      IF (INIT.EQ.'R' .OR. INIT.EQ.'r') THEN
         IF (NINT(TRIG(N)).NE.N .OR. NINT(TRIG(2*N)).NE.N) THEN
            IERROR = 5
            RETURN
         END IF
      END IF
      CALL C06FPZ(N,NQ,Q)
      IF (INIT.EQ.'I' .OR. INIT.EQ.'i') THEN
         CALL C06FPY(N,NQ,Q,TRIG(1),TRIG(N+1))
      END IF
      NCHECK = N
      RETURN
      END
