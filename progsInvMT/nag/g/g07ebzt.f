      SUBROUTINE G07EBZ(LOWER,ALPHA,N,M,IND,WRK,LWRK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA
      INTEGER           IND, LWRK, M, N
      LOGICAL           LOWER
C     .. Array Arguments ..
      DOUBLE PRECISION  WRK(LWRK)
C     .. Local Scalars ..
      DOUBLE PRECISION  CC, PH0, PH1, UMU, USIGMA, X0, X1, XI0, XI1, XM,
     *                  XN, Z
      INTEGER           I0, I1, IF2, NM, NMMAX
      LOGICAL           EXACT, PH1SET
C     .. External Functions ..
      DOUBLE PRECISION  G01EAF, G01FAF
      EXTERNAL          G01EAF, G01FAF
C     .. External Subroutines ..
      EXTERNAL          G08AJF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, NINT, SQRT
C     .. Executable Statements ..
      PH1SET = .FALSE.
      NMMAX = MAX(N,M)
      XN = DBLE(N)
      XM = DBLE(M)
      NM = N*M
      EXACT = .FALSE.
      IF ((N+M).LE.40 .OR. NMMAX.LE.30) EXACT = .TRUE.
C
C     Find an initial estimate based on the normal approximation.
C
      UMU = DBLE(NM)
      USIGMA = SQRT(XN*XM*(XN+XM+1.0D0)/3.0D0)
      IF2 = 0
      IF (LOWER) THEN
         CC = 1.0D0
         Z = G01FAF('Lower',ALPHA,IF2)
         XI0 = (Z*USIGMA+UMU-CC)/2.0D0
         I0 = NINT(XI0)
         IF (I0.LT.-1) I0 = -1
      ELSE
         CC = -1.0D0
         IF2 = 0
         Z = G01FAF('Upper',ALPHA,IF2)
         XI0 = (Z*USIGMA+UMU-CC)/2.0D0
         I0 = NINT(XI0)
         IF (I0.GT.NM+1) I0 = NM + 1
      END IF
C
C     Compute the tail probability PH0 corresponding to I0.
C
   20 CONTINUE
      IF (I0.LT.0) THEN
         PH0 = 0.0D0
      ELSE IF (I0.GT.NM) THEN
         PH0 = 0.0D0
      ELSE
         XI0 = DBLE(I0)
         IF (EXACT) THEN
C
C           Exact significance level for small samples.
C
            IF (LOWER) THEN
               CALL G08AJF(N,M,'Lower',XI0,PH0,WRK,LWRK,IF2)
            ELSE
               CALL G08AJF(N,M,'Upper',XI0,PH0,WRK,LWRK,IF2)
            END IF
         ELSE
            IF2 = 0
            X0 = (2.0D0*XI0-UMU+CC)/USIGMA
            IF (LOWER) THEN
               PH0 = G01EAF('Lower',X0,IF2)
            ELSE
               PH0 = G01EAF('Upper',X0,IF2)
            END IF
         END IF
      END IF
   40 CONTINUE
C
C     Check whether PH0 satisfies the required probability ALPHA
C
      IF (PH0.GT.ALPHA) THEN
C
C        If not set I1 to I0 and shift I0 by one. Go back and
C        recompute PH0.
C
         I1 = I0
         PH1 = PH0
         PH1SET = .TRUE.
         IF (LOWER) THEN
            I0 = I0 - 1
         ELSE
            I0 = I0 + 1
         END IF
         GO TO 20
      ELSE IF ( .NOT. PH1SET) THEN
C
C        PH1 must still be calculated.
C
         IF (LOWER) THEN
            I1 = I0 + 1
            XI1 = DBLE(I1)
            IF (EXACT) THEN
               CALL G08AJF(N,M,'LOWER',XI1,PH1,WRK,LWRK,IF2)
            ELSE
               IF2 = 0
               X1 = (2.0D0*XI1-UMU+CC)/USIGMA
               PH1 = G01EAF('Lower',X1,IF2)
            END IF
         ELSE
            I1 = I0 - 1
            XI1 = DBLE(I1)
            IF (EXACT) THEN
               CALL G08AJF(N,M,'UPPER',XI1,PH1,WRK,LWRK,IF2)
            ELSE
               IF2 = 0
               X1 = (2.0D0*XI1-UMU+CC)/USIGMA
               PH1 = G01EAF('Upper',X1,IF2)
            END IF
         END IF
      END IF
C
C     Now check PH1. If PH1 is .le. ALPHA then need to shift I0 and I1
C     by one, which means PH1 needs to be recalculated.
C
      IF (PH1.LE.ALPHA) THEN
         I0 = I1
         PH0 = PH1
         PH1SET = .FALSE.
         GO TO 40
      END IF
C
C     Finished - set IND to I0.
C
      IND = I0
      IF (IND.LT.0) THEN
         IND = 0
      ELSE IF (IND.GT.NM) THEN
         IND = NM
      END IF
C
      RETURN
      END
