      SUBROUTINE F08MET(N,Q,E,TAU,SUP)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Purpose
C     =======
C
C     DLASQ4 estimates TAU, the smallest eigenvalue of B'B
C     where B is a upper bidiagonal matrix. Q and E contains the
C     squared values of the diagonal and the superdiagonal respectively.
C     This routine improves the input value of SUP which is an upper
C     bound for the smallest eigenvalue for this matrix.
C
C     Arguments
C     =========
C
C  N       (input) INTEGER
C          On entry, N specifies the number of rows and columns
C          in the matrix. N must be at least 0.
C
C  Q       (input) DOUBLE PRECISION array, dimension (N)
C          Q array
C
C  E       (input) DOUBLE PRECISION array, dimension (N-1)
C          E array
C
C  TAU     (output) DOUBLE PRECISION
C          Estimate of the shift
C
C  SUP     (input/output) DOUBLE PRECISION
C          Upper bound for the smallest singular value
C
C  -- LAPACK routine (version 2.0) (adapted for NAG library) --
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C     September 30, 1994
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
      DOUBLE PRECISION  BIS, BIS1
      PARAMETER         (BIS=0.9999D+0,BIS1=0.7D+0)
      INTEGER           IFLMAX
      PARAMETER         (IFLMAX=5)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SUP, TAU
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  E(*), Q(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  D, DM, XINF
      INTEGER           I, IFL
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
      IFL = 1
      SUP = MIN(SUP,Q(1),Q(2),Q(3),Q(N),Q(N-1),Q(N-2))
      TAU = SUP*BIS
      XINF = ZERO
   20 CONTINUE
      IF (IFL.EQ.IFLMAX) THEN
         TAU = XINF
         RETURN
      END IF
      D = Q(1) - TAU
      DM = D
      DO 40 I = 1, N - 2
         D = (D/(D+E(I)))*Q(I+1) - TAU
         IF (DM.GT.D) DM = D
         IF (D.LT.ZERO) THEN
            SUP = TAU
            TAU = MAX(SUP*BIS1**IFL,D+TAU)
            IFL = IFL + 1
            GO TO 20
         END IF
   40 CONTINUE
      D = (D/(D+E(N-1)))*Q(N) - TAU
      IF (DM.GT.D) DM = D
      IF (D.LT.ZERO) THEN
         SUP = TAU
         XINF = MAX(XINF,D+TAU)
         IF (SUP*BIS1**IFL.LE.XINF) THEN
            TAU = XINF
         ELSE
            TAU = SUP*BIS1**IFL
            IFL = IFL + 1
            GO TO 20
         END IF
      ELSE
         SUP = MIN(SUP,DM+TAU)
      END IF
      RETURN
C
C     End of F08MET (DLASQ4)
C
      END
