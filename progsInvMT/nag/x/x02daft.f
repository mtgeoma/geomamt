      LOGICAL FUNCTION X02DAF(X)
C     MARK 8 RELEASE. NAG COPYRIGHT 1980.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     RETURNS .FALSE. IF THE SYSTEM SETS UNDERFLOWING QUANTITIES
C     TO ZERO, WITHOUT ANY ERROR INDICATION OR UNDESIRABLE WARNING
C     OR SYSTEM OVERHEAD.
C     RETURNS .TRUE. OTHERWISE, IN WHICH CASE CERTAIN LIBRARY
C     ROUTINES WILL TAKE SPECIAL PRECAUTIONS TO AVOID UNDERFLOW
C     (USUALLY AT SOME COST IN EFFICIENCY).
C
C     X IS A DUMMY ARGUMENT
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION        X
C     .. Executable Statements ..
      X02DAF = .FALSE.
      RETURN
      END