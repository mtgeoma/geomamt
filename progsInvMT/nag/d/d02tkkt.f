      SUBROUTINE D02TKK(UHIGH,HI,DMZ,NEQ,KCOL,COEF)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C**********************************************************************
C
C   Purpose:
C      Determine highest order (piecewise constant) derivatives
C      of the current collocation solution
C
C   Arguments:
C      HI     - the stepsize, HI = XI(I+1) - XI(I), for the
C               interval being treated
C      DMZ    - vector of mj-th derivative of the solution at each
C               collocation point of this interval
C      UHIGH  - the array of highest order (piecewise constant)
C               derivatives of the approximate solution on
C               (xi(i),xi(i+1)), viz,
C                             (k+mj-1)
C                 UHIGH(j) = u   (x)    on (xi(i),xi(i+1));
C                             j
C               each entry is a linear combination of the corresponding
c               entries in DMZ
C      NEQ    - number of ODEs
C      KCOL   - number of collocation points per interval
C      COEF   - rk-basis coefficients
C
C   Author:
C      R.W. Brankin, NAG Ltd., August 1994
C      (modified version of COLNEW routine HORDER)
C
C**********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HI
      INTEGER           KCOL, NEQ
C     .. Array Arguments ..
      DOUBLE PRECISION  COEF(KCOL,KCOL), DMZ(NEQ,KCOL), UHIGH(NEQ)
C     .. Local Scalars ..
      DOUBLE PRECISION  DN, FACT
      INTEGER           ID, J, KIN
C     .. Executable Statements ..
C
      DN = 1.D0/HI**(KCOL-1)
C
C Loop over the NEQ solution components
C
      DO 20 ID = 1, NEQ
         UHIGH(ID) = 0.D0
   20 CONTINUE
      KIN = 1
      DO 60 J = 1, KCOL
         FACT = DN*COEF(1,J)
         DO 40 ID = 1, NEQ
            UHIGH(ID) = UHIGH(ID) + FACT*DMZ(ID,J)
   40    CONTINUE
c         KIN = KIN + Kcol
   60 CONTINUE
      RETURN
      END
