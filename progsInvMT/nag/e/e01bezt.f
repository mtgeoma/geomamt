      SUBROUTINE E01BEZ(N,X,F,D,INCFD,IERR)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ------------------------------------------------------------------
C
C     Derived from PCHIP routine
C          PCHIM : Piecewise Cubic Hermite Interpolation to
C                  Monotone Data.
C
C     Sets derivatives needed to determine a monotone piecewise cubic
C     Hermite interpolant to the data given in X and F.
C
C     Default boundary conditions are provided which are compatible
C     with monotonicity.
C
C     If the data are only piecewise monotonic, the interpolant will
C     have an extremum at each point where monotonicity switches direc-
C     tion.
C
C     To facilitate two-dimensional applications, includes an increment
C     between successive values of the F- and D-arrays.
C
C     The resulting piecewise cubic Hermite function may be evaluated
C     by E01BFF or E01BGF.
C
C     ------------------------------------------------------------------
C
C     Parameters:
C
C     N -- (input) number of data points. N.ge.2.
C           If N=2, simply does linear interpolation.
C
C     X -- (input) real array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .lt. X(I),  I = 2(1)N.
C
C     F -- (input) real array of dependent variable values to be inter-
C           polated.  F(1+(I-1)*INCFD) is value corresponding to X(I).
C           E01BEF is designed for monotonic data, but it will work for
C           any F-array.  It will force extrema at points where mono-
C           tonicity switches direction.
C
C     D -- (output) real array of derivative values at the data points.
C           If the data are monotonic, these values will determine a
C           a monotone cubic Hermite function.
C           The value corresponding to X(I) is stored in
C                D(1+(I-1)*INCFD),  I=1(1)N.
C           No other entries in D are changed.
C
C     INCFD -- (input) increment between successive values in F and D.
C           This argument is provided primarily for 2-d applications.
C           INCFD.ge.1.
C
C     ------------------------------------------------------------------
C
C     References:  1. F.N.Fritsch and R.E.Carlson, 'Monotone Piecewise
C                     Cubic Interpolation', SIAM J. Numer. Anal. 17, 2
C                     (April 1980), 238-246.
C                  2. F.N.Fritsch and J.Butland, 'A Method for
C                     Constructing Local Monotone Piecewise Cubic
C                     Interpolants',
C                     LLNL Preprint UCRL-87559 (April 1982).
C
C     ------------------------------------------------------------------
C
C     Programmed by:  Fred N. Fritsch,  FTS 532-4275,
C                     Mathematics and Statistics Division,
C                     Lawrence Livermore National Laboratory.
C
C     Change record:
C     82-02-01   1. Introduced E01BEY to reduce possible over/under-
C                   flow problems.
C                2. Rearranged derivative formula for same reason.
C     82-06-02   1. Modified end conditions to be continuous functions
C                   of data when monotonicity switches in next interval.
C                2. Modified formulas so end conditions are less prone
C                   to over/underflow problems.
C     82-08-03      Minor cosmetic changes for release 1.
C
C     ------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, THREE
      PARAMETER         (ZERO=0.0D0,THREE=3.0D0)
C     .. Scalar Arguments ..
      INTEGER           IERR, INCFD, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(INCFD,N), F(INCFD,N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  DEL1, DEL2, DMAX, DMIN, DRAT1, DRAT2, DSAVE, H1,
     *                  H2, HSUM, HSUMT3, TMP, W1, W2
      INTEGER           I, NLESS1
C     .. External Functions ..
      DOUBLE PRECISION  E01BEY
      EXTERNAL          E01BEY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Executable Statements ..
C     Assume that arguments have been checked by E01BEF.
C
      IERR = 0
      NLESS1 = N - 1
      H1 = X(2) - X(1)
      DEL1 = (F(1,2)-F(1,1))/H1
      DSAVE = DEL1
C     Special case N=2 -- use linear interpolation.
      IF (NLESS1.LE.1) THEN
         D(1,1) = DEL1
         D(1,N) = DEL1
         GO TO 40
      END IF
C     Normal case  (N .ge. 3).
      H2 = X(3) - X(2)
      DEL2 = (F(1,3)-F(1,2))/H2
C     Set D(1) via non-centered three-point formula, adjusted to be
C     shape-preserving.
      HSUM = H1 + H2
      W1 = (H1+HSUM)/HSUM
      W2 = -H1/HSUM
      D(1,1) = W1*DEL1 + W2*DEL2
      IF (E01BEY(D(1,1),DEL1).LE.ZERO) THEN
         D(1,1) = ZERO
      ELSE IF (E01BEY(DEL1,DEL2).LT.ZERO) THEN
C        need do this check only if monotonicity switches.
         DMAX = THREE*DEL1
         IF (ABS(D(1,1)).GT.ABS(DMAX)) D(1,1) = DMAX
      END IF
C     Loop through interior points.
      DO 20 I = 2, NLESS1
         IF (I.NE.2) THEN
            H1 = H2
            H2 = X(I+1) - X(I)
            HSUM = H1 + H2
            DEL1 = DEL2
            DEL2 = (F(1,I+1)-F(1,I))/H2
         END IF
C        Set D(I)=0 unless data are strictly monotonic.
         D(1,I) = ZERO
         TMP = E01BEY(DEL1,DEL2)
         IF (TMP.LT.ZERO) THEN
            IERR = IERR - 1
            DSAVE = DEL2
         ELSE IF (TMP.EQ.ZERO) THEN
C           Count number of changes in direction of monotonicity.
            IF (DEL2.NE.ZERO) THEN
               IF (E01BEY(DSAVE,DEL2).LT.ZERO) IERR = IERR - 1
               DSAVE = DEL2
            END IF
         ELSE
C           Use Brodlie modification of Butland formula.
            HSUMT3 = HSUM + HSUM + HSUM
            W1 = (HSUM+H1)/HSUMT3
            W2 = (HSUM+H2)/HSUMT3
            DMAX = MAX(ABS(DEL1),ABS(DEL2))
            DMIN = MIN(ABS(DEL1),ABS(DEL2))
            DRAT1 = DEL1/DMAX
            DRAT2 = DEL2/DMAX
            D(1,I) = DMIN/(W1*DRAT1+W2*DRAT2)
         END IF
   20 CONTINUE
C
C     Set D(N) via non-centered three-point formula, adjusted to be
C     shape-preserving.
      W1 = -H2/HSUM
      W2 = (H2+HSUM)/HSUM
      D(1,N) = W1*DEL1 + W2*DEL2
      IF (E01BEY(D(1,N),DEL2).LE.ZERO) THEN
         D(1,N) = ZERO
      ELSE IF (E01BEY(DEL1,DEL2).LT.ZERO) THEN
C        Need do this check only if monotonicity switches.
         DMAX = THREE*DEL2
         IF (ABS(D(1,N)).GT.ABS(DMAX)) D(1,N) = DMAX
      END IF
   40 CONTINUE
C
C     Normal return.
C
C     Negative IERR signifies -IERR changes in monotonicity direction.
C
      RETURN
      END
