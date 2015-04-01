      DOUBLE PRECISION FUNCTION E01BHX(X1,X2,F1,F2,D1,D2,A,B,IERR)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ------------------------------------------------------------------
C
C     Derived from PCHIP routine
C          CHFIV:  Cubic Hermite Function Integral Evaluator.
C
C     Called by  E01BHF  to evaluate the integral of a single cubic (in
C     Hermite form) over an arbitrary interval (A,B).
C
C     ------------------------------------------------------------------
C
C     Parameters:
C
C     VALUE -- (output) value of the requested integral.
C
C     X1,X2 -- (input) endpoints if interval of definition of cubic.
C              (Must be distinct.  Error return if not.)
C
C     F1,F2 -- (input) function values at the ends of the interval.
C
C     D1,D2 -- (input) derivative values at the ends of the interval.
C
C     A,B --   (input) endpoints of interval of integration.
C
C     IERR --  (output) error flag.
C              IERR = 0  Normal return.
C              IERR = 1  if X1 = X2.
C              (Value has not been set if IERR = 1.)
C
C     ------------------------------------------------------------------
C
C     Programmed by:  Fred N. Fritsch,  FTS 532-4275,
C                     Mathematics and Statistics Division,
C                     Lawrence Livermore National Laboratory.
C
C     ------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION                 HALF, TWO, THREE, FOUR, SIX
      PARAMETER                        (HALF=0.5D0,TWO=2.0D0,
     *                                 THREE=3.0D0,FOUR=4.0D0,SIX=6.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, B, D1, D2, F1, F2, X1, X2
      INTEGER                          IERR
C     .. Local Scalars ..
      DOUBLE PRECISION                 DTERM, FTERM, H, PHIA1, PHIA2,
     *                                 PHIB1, PHIB2, PSIA1, PSIA2,
     *                                 PSIB1, PSIB2, TA1, TA2, TB1, TB2,
     *                                 UA1, UA2, UB1, UB2
C     .. Executable Statements ..
C     Validity check input.
      IF (X1.EQ.X2) THEN
C
C        Error return.
         IERR = 1
      ELSE
         IERR = 0
C        Compute integral.
         H = X2 - X1
         TA1 = (A-X1)/H
         TA2 = (X2-A)/H
         TB1 = (B-X1)/H
         TB2 = (X2-B)/H
C
         UA1 = TA1**3
         PHIA1 = UA1*(TWO-TA1)
         PSIA1 = UA1*(THREE*TA1-FOUR)
         UA2 = TA2**3
         PHIA2 = UA2*(TWO-TA2)
         PSIA2 = -UA2*(THREE*TA2-FOUR)
C
         UB1 = TB1**3
         PHIB1 = UB1*(TWO-TB1)
         PSIB1 = UB1*(THREE*TB1-FOUR)
         UB2 = TB2**3
         PHIB2 = UB2*(TWO-TB2)
         PSIB2 = -UB2*(THREE*TB2-FOUR)
C
         FTERM = F1*(PHIA2-PHIB2) + F2*(PHIB1-PHIA1)
         DTERM = (D1*(PSIA2-PSIB2)+D2*(PSIB1-PSIA1))*(H/SIX)
C        Return value.
         E01BHX = (HALF*H)*(FTERM+DTERM)
      END IF
      RETURN
      END
