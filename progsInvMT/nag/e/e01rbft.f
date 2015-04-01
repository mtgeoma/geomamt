      SUBROUTINE E01RBF(M,A,U,X,F,IFAIL)
C     MARK 9 RELEASE. NAG COPYRIGHT 1981
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     E01RBF EVALUATES CONTINUED FRACTIONS OF THE FORM
C     PRODUCED BY E01RAF
C
C     ON ENTRY
C     X  -REAL-  POINT AT WHICH THIELE INTERPOLANT IS TO BE
C     EVALUATED
C     U  -REAL ARRAY-  ABSCISSA POINTS (REARRANGED)
C     A  -REAL ARRAY-  THIELE COEFFICIENTS
C     M  -INTEGER-  NUMBER OF COEFFICIENTS
C     ON EXIT
C     F  -REAL-  CALCULATED VALUE OF THE THIELE INTERPOLANT
C     IFAIL  -INTEGER-  ERROR INDICATOR
C
C     EPS IS A MACHINE DEPENDENT VARIABLE.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E01RBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  F, X
      INTEGER           IFAIL, M
C     .. Array Arguments ..
      DOUBLE PRECISION  A(M), U(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, W, W1, W2, W3, W4, W5, W6
      INTEGER           IFL, J, MM1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      EPS = SQRT(X02AJF())
      IFL = IFAIL
      IFAIL = 0
      IF (M.NE.1) GO TO 20
      F = A(1)
      RETURN
   20 MM1 = M - 1
      W = 0.0D0
      W1 = 1.0D0
      W2 = 1.0D0
      W3 = A(1)
      DO 40 J = 1, MM1
         W6 = A(J+1)*(X-U(J))
         W4 = W1 + W6*W
         W5 = W3 + W6*W2
         W = W1
         W1 = W4
         W2 = W3
         W3 = W5
   40 CONTINUE
C     TEST FOR A ZERO DENOMINATOR.
      IF (ABS(W4).GE.EPS) GO TO 60
      IFAIL = P01ABF(IFL,1,SRNAME,0,P01REC)
      RETURN
   60 F = W5/W4
      RETURN
      END
