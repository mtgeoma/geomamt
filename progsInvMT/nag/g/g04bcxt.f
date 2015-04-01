      SUBROUTINE G04BCX(N,NREP,NROW,NCOL,NT,IT,C,LDC,IREP,EF,NTDF,TMEAN,
     *                  ACC,WK,IWARN)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Computes treatment effects for non-orthogonal row-column designs
C     and the generalised inverse for the reduced treatment matrix
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ACC
      INTEGER           IWARN, LDC, N, NCOL, NREP, NROW, NT, NTDF
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LDC,NT), EF(NT), TMEAN(NT), WK(3*NT)
      INTEGER           IREP(NT), IT(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  RMEAN, VAR
      INTEGER           I, ICOL, IFAULT, J
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, F02FAF, G04BCY
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
      IWARN = 0
      DO 40 I = 1, NT
         DO 20 J = 1, NT
            C(J,I) = 0.0D0
   20    CONTINUE
   40 CONTINUE
C
C      Compute the NN' matrices
C
      CALL G04BCY(NREP,NROW*NCOL,NT,IT,C,LDC,-1.0D0)
      CALL G04BCY(NREP*NROW,NCOL,NT,IT,C,LDC,1.0D0)
      DO 60 I = 1, NREP
         CALL G04BCY(-NCOL,NROW,NT,IT((I-1)*NROW*NCOL+1),C,LDC,1.0D0)
   60 CONTINUE
C
C     Compute A matrix
C
      DO 80 I = 1, NT
         C(I,I) = IREP(I) + C(I,I)
   80 CONTINUE
C
C     Find efficiency factors and parameter estimates
C
C     Compute eigenvalues of A
C
      IFAULT = 1
      CALL F02FAF('V','L',NT,C,LDC,EF,WK,3*NT,IFAULT)
      IF (IFAULT.NE.0) THEN
         IWARN = -4
         GO TO 240
      END IF
C
C     Check for zero values and compute q*U*INV(E)
C
      NTDF = 0
      DO 100 I = 1, NT
         IF (EF(I).GT.ACC) THEN
            NTDF = NTDF + 1
            WK(I) = DDOT(NT,C(1,I),1,TMEAN,1)/EF(I)
         ELSE
            EF(I) = 0.0D0
         END IF
  100 CONTINUE
      IF (NTDF.EQ.0) THEN
         IWARN = -3
         GO TO 240
      END IF
      ICOL = NT - NTDF + 1
      CALL DGEMV('N',NT,NTDF,1.0D0,C(1,ICOL),LDC,WK(ICOL),1,0.0D0,TMEAN,
     *           1)
C
C     Compute generalised inverse
C
      DO 140 I = NT, 1, -1
         DO 120 J = ICOL, NT
            WK(J) = C(I,J)/EF(J)
  120    CONTINUE
         CALL DGEMV('N',I,NTDF,1.0D0,C(1,ICOL),LDC,WK(ICOL),1,0.0D0,
     *              WK(NT+1),1)
         CALL DCOPY(I,WK(NT+1),1,C(I,1),LDC)
  140 CONTINUE
      DO 160 I = 1, NT - 1
         CALL DCOPY(NT-I,C(I+1,I),1,C(I,I+1),LDC)
  160 CONTINUE
C
C     Compute se of differences in means
C
      DO 200 I = 1, NT - 1
         DO 180 J = I + 1, NT
            VAR = C(I,I) + C(J,J) - 2.0D0*C(I,J)
            IF (VAR.GT.0.0D0) THEN
               C(J,I) = SQRT(VAR)
            ELSE
               C(J,I) = 0.0D0
               IWARN = -1
            END IF
  180    CONTINUE
  200 CONTINUE
C
C     Scale efficiency factors
C
      RMEAN = DBLE(N)/DBLE(NT)
      DO 220 I = ICOL, NT
         EF(I) = EF(I)/RMEAN
  220 CONTINUE
      IF (NTDF.LT.NT-1) IWARN = -2
  240 CONTINUE
      RETURN
      END
