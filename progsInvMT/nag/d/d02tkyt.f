      SUBROUTINE D02TKY(S,COEF,KCOL,MAXORD,RKB,GETDMB,DM)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C**********************************************************************
C
C   Purpose:
C      Evaluate mesh independent runge-kutta basis for given S
C
C   Variables:
C     S      - argument, i.e. the relative position for which
C              the basis is to be evaluated ( 0 <= S <= 1 ).
C     COEF   - precomputed derivatives of the basis
C     KCOL   - number of collocation points per subinterval
C     MAXORD - maximal order of the differential equation
C     RKB    - the runge-kutta basis (0-th to (m-1)-th derivatives )
C     GETDMB - determine if basis for maximal order to be computed
C     DM     - basis elements for m-th derivative
C
C
C    Author:
C       R.W. Brankin, NAG Ltd., August 1994
C       (modified version of the COLNEW routine RKBAS)
C
C**********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION  S
      INTEGER           KCOL, MAXORD
      LOGICAL           GETDMB
C     .. Array Arguments ..
      DOUBLE PRECISION  COEF(KCOL,KCOL), DM(KCOL), RKB(7,MAXORD)
C     .. Local Scalars ..
      DOUBLE PRECISION  P
      INTEGER           I, J, KPM1, L, LB
C     .. Local Arrays ..
      DOUBLE PRECISION  T(10)
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
      IF (KCOL.EQ.1) THEN
         RKB(1,1) = 1.0D0
         DM(1) = 1.0D0
      ELSE
         KPM1 = KCOL + MAXORD - 1
         DO 20 I = 1, KPM1
            T(I) = S/DBLE(I)
   20    CONTINUE
         DO 80 L = 1, MAXORD
            LB = KCOL + L + 1
            DO 60 I = 1, KCOL
               P = COEF(1,I)
               DO 40 J = 2, KCOL
                  P = P*T(LB-J) + COEF(J,I)
   40          CONTINUE
               RKB(I,L) = P
   60       CONTINUE
   80    CONTINUE
         IF (GETDMB) THEN
            DO 120 I = 1, KCOL
               P = COEF(1,I)
               DO 100 J = 2, KCOL
                  P = P*T(KCOL+1-J) + COEF(J,I)
  100          CONTINUE
               DM(I) = P
  120       CONTINUE
         END IF
      END IF
      RETURN
      END
