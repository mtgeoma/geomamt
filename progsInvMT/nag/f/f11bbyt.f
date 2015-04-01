      SUBROUTINE F11BBY(NEXT,IREVCM,NORM,ANORM,N,WORK,U,V)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-----------------------------------------------------------------------
C
C     F11BBY - Compute the norm of the matrix A
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER           ITMAX
      PARAMETER         (ITMAX=5)
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ANORM
      INTEGER           IREVCM, N, NORM
      LOGICAL           NEXT
C     .. Array Arguments ..
      DOUBLE PRECISION  U(N), V(N), WORK(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALTSGN, ESTOLD
      INTEGER           I, ITER, J, JLAST, JUMP, KASE
      LOGICAL           DONE, VARSIG
C     .. External Functions ..
      DOUBLE PRECISION  DASUM
      INTEGER           IDAMAX
      EXTERNAL          DASUM, IDAMAX
C     .. External Subroutines ..
      EXTERNAL          F06FBF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, MOD, NINT, SIGN
C     .. Save statement ..
      SAVE              ITER, J, JUMP, KASE
C     .. Executable Statements ..
C
C     First call to this routine
C
      IF (NEXT) THEN
         NEXT = .FALSE.
         ANORM = ZERO
         IF (MOD(NORM,3).LE.0) THEN
            IREVCM = 1
         ELSE
            IREVCM = -1
         END IF
         CALL F06FBF(N,ONE/DBLE(N),U,1)
         KASE = 1
         JUMP = 1
         ANORM = ZERO
C
C     Subsequent calls to this routine
C
      ELSE
         DONE = .FALSE.
C
C        First iteration (JUMP = 1): u overwritten by A*u
C
         IF (JUMP.LE.1) THEN
            ANORM = DASUM(N,V,1)
            IF (N.EQ.1) THEN
               KASE = 0
            ELSE
               DO 20 I = 1, N
                  U(I) = SIGN(ONE,V(I))
                  WORK(I) = NINT(U(I))
   20          CONTINUE
               KASE = 2
               JUMP = 2
            END IF
C
C        First iteration (JUMP = 2): u overwritten by A'*u
C
         ELSE IF (JUMP.EQ.2) THEN
            ITER = 2
            J = IDAMAX(N,V,1)
            CALL F06FBF(N,ZERO,U,1)
            U(J) = ONE
            KASE = 1
            JUMP = 3
C
C        Subsequent iteration (JUMP = 3): u overwritten by A*u
C
         ELSE IF (JUMP.EQ.3) THEN
            ESTOLD = ANORM
            DO 40 I = 1, N
               VARSIG = NINT(SIGN(ONE,V(I))) .NE. NINT(WORK(I))
               IF (VARSIG) GO TO 60
   40       CONTINUE
   60       CONTINUE
            ANORM = DASUM(N,V,1)
C
            IF (( .NOT. VARSIG) .OR. (ANORM.LE.ESTOLD)) THEN
               DONE = .TRUE.
            ELSE
               DO 80 I = 1, N
                  U(I) = SIGN(ONE,V(I))
                  WORK(I) = NINT(U(I))
   80          CONTINUE
               KASE = 2
               JUMP = 4
            END IF
C
C        Subsequent iteration (JUMP = 4): u overwritten by A'*u
C
         ELSE IF (JUMP.EQ.4) THEN
            JLAST = J
            J = IDAMAX(N,V,1)
            IF ((ABS(V(J)).NE.ABS(V(JLAST))) .AND. (ITER.LT.ITMAX)) THEN
               ITER = ITER + 1
               CALL F06FBF(N,ZERO,U,1)
               U(J) = ONE
               KASE = 1
               JUMP = 3
            ELSE
               DONE = .TRUE.
            END IF
C
C        Subsequent iteration (JUMP = 5): u overwritten by A'*u
C
         ELSE
            ANORM = MAX(ANORM,TWO*DASUM(N,V,1)/DBLE(3*N))
            KASE = 0
         END IF
C
C        Iteration complete: final stage
C
         IF (DONE) THEN
            ALTSGN = 1
            DO 100 I = 1, N
               U(I) = ALTSGN*(1+DBLE(I-1)/DBLE(N-1))
               ALTSGN = -ALTSGN
  100       CONTINUE
            KASE = 1
            JUMP = 5
         END IF
C
C        Completion
C
         IF (KASE.NE.0) THEN
            IREVCM = -IREVCM
         ELSE
            NEXT = .TRUE.
         END IF
C
      END IF
C
C     End of subroutine F11BBY
C
      RETURN
      END
