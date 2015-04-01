      SUBROUTINE G05EGF(E,A,NA,B,NB,R,NR,VAR,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 14 REVISED. IER-747 (DEC 1989).
C
C     G05EGF SETS UP A REFERENCE VECTOR FOR AN AUTOREGRESSIVE
C     MOVING-AVERAGE TIME-SERIES MODEL WITH NORMALLY
C     DISTRIBUTED ERRORS, SO THAT G05EWF MAY BE USED TO
C     GENERATE SUCCESSIVE TERMS. IT ALSO INITIALISES
C     THE  SERIES TO A STATIONARY POSITION.
C
C     THE METHOD IS FROM A PAPER OF G.TUNNICLIFFE WILSON IN THE
C     JOURNAL FOR STATISTICAL COMPUTATION AND SIMULATION (C.1978).
C     THE REFERENCE VECTOR CONTAINS THE FOLLOWING VALUES
C     1) NA, THE AUTOREGRESSIVE ORDER.
C     2) NB, THE MOVING AVERAGE ORDER.
C     3) A POINTER TO THE LAST VALUE OF THE PURE AUTOREGRESSIVE
C        SERIES.
C     4) E, THE MEAN OF THE ERROR TERM.
C     5) A, THE AUTOREGRESSIVE COEFFICIENTS.
C     NA+5) B, THE MOVING-AVERAGE COEFFICIENTS.
C     NA+NB+5) THE LAST MAX(NA,NB) VALUES OF THE PURE AUTOREGRESSIVE
C              SERIES IN A CIRCULAR BUFFER.
C
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G05EGF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  E, VAR
      INTEGER           IFAIL, NA, NB, NR
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*), B(NB), R(NR)
C     .. Local Scalars ..
      DOUBLE PRECISION  ERR, HALF, ONE, TWO, X, Y, Z, ZERO
      INTEGER           I, IERR, II, J, JJ, K, L, MMAX
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G05DDF, X02AJF
      INTEGER           P01ABF
      EXTERNAL          G05DDF, X02AJF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, DBLE, SQRT
C     .. Data statements ..
      DATA              ZERO/0.0D0/, HALF/0.5D0/, ONE/1.0D0/, TWO/2.0D0/
C     .. Executable Statements ..
C
C     CHECK ARGUMENTS FOR SENSIBLE VALUES AND CONSISTENCY
C
      IERR = 1
      IF (NA.LT.0) GO TO 280
      IERR = 2
      IF (NB.LT.1) GO TO 280
      IERR = 3
      MMAX = NA + NB + 4 + MAX(NA,NB)
      IF (NR.LT.MMAX) GO TO 280
      ERR = TWO*DBLE(NA+1)*X02AJF()
      IERR = 4
C
C     DECOMPOSE THE AUTOREGRESSION INTO ITS LOWER ORDER ONES, USING
C     THE REFERENCE VECTOR AS WORKSPACE.
C
      VAR = ONE
      L = MMAX
      IF (NB.GT.NA) L = L - NB + NA
      DO 20 I = 1, NA
         R(I) = A(I)
   20 CONTINUE
      DO 80 II = 1, NA
         I = NA + 1 - II
         Z = ONE - R(I)**2
         VAR = VAR*Z
         IF (VAR.LT.ERR) GO TO 280
         JJ = (I-1)/2
         IF (I-1.GT.2*JJ) R(JJ+1) = (ONE+R(I))*R(JJ+1)/Z
         IF (JJ.LT.1) GO TO 60
         DO 40 J = 1, JJ
            X = R(J)
            K = I - J
            Y = R(K)
            R(J) = (X+R(I)*Y)/Z
            R(K) = (Y+R(I)*X)/Z
   40    CONTINUE
   60    R(L) = G05DDF(ZERO,ONE/SQRT(VAR))
         L = L - 1
   80 CONTINUE
      IF (NA.EQ.0) X = G05DDF(ZERO,ONE)
C
C     APPLY THE AUTOREGRESSIONS IN TURN IN ORDER TO BUILD UP THE
C     VECTOR.
C
      L = L + 1
      II = NA - 1
      IF (II.LT.1) GO TO 160
      DO 140 I = 1, II
         L = L + 1
         DO 100 K = 1, I
            J = L - K
            R(L) = R(L) + R(J)*R(K)
  100    CONTINUE
         JJ = I/2
         IF (I.GT.2*JJ) R(JJ+1) = (ONE-R(I+1))*R(JJ+1)
         IF (JJ.LT.1) GO TO 140
         DO 120 J = 1, JJ
            X = R(J)
            K = I + 1 - J
            Y = R(K)
            R(J) = X - R(I+1)*Y
            R(K) = Y - R(I+1)*X
  120    CONTINUE
  140 CONTINUE
  160 IF (L.GE.MMAX) GO TO 220
      L = L + 1
      DO 200 I = L, MMAX
         X = G05DDF(ZERO,ONE)
         J = I
         DO 180 K = 1, NA
            J = J - 1
            X = X + R(J)*A(K)
  180    CONTINUE
         R(I) = X
  200 CONTINUE
C
C     COPY THE DIMENSIONS AND INPUT ARRAYS.
C
  220 R(1) = DBLE(NA) + HALF
      R(2) = DBLE(NB) + HALF
      R(3) = DBLE(MMAX) + HALF
      R(4) = E
      DO 240 I = 1, NA
         R(I+4) = A(I)
  240 CONTINUE
      DO 260 I = 1, NB
         II = NA + 4 + I
         R(II) = B(I)
  260 CONTINUE
      IFAIL = 0
      RETURN
  280 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
      END
