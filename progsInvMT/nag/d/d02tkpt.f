      SUBROUTINE D02TKP(H,B,GI,MSTAR,VI,M,MAXORD,NEQ,KCOL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C**********************************************************************
C
C   Purpose:
C
C      Construct collocation matrix rows for an interior mesh interval.
C
C   Arguments:
C
C      H      - local stepsize
C      B      - rk-basis coefficients
C      GI     - sub-block of the collocation matrix in
C               which the equations are to be formed
C      MSTAR  - number of unknowns in z(u(x))
C      VI     - sub-block of noncondensed collocation equations,
C               right-hand side part
C      M      - orders of ODEs
C      MAXORD - maximal order of ODEs
C      NEQ    - number of ODEs
C      KCOL   - number of collocation points
C
C   Author:
C      R.W. Brankin, NAG Ltd., August 1994
C      (modified version of part of COLNEW routine GBLOCK)
C
C**********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H
      INTEGER           KCOL, MAXORD, MSTAR, NEQ
C     .. Array Arguments ..
      DOUBLE PRECISION  B(7,4), GI(MSTAR,2*MSTAR), VI(NEQ,KCOL,MSTAR)
      INTEGER           M(NEQ)
C     .. Local Scalars ..
      DOUBLE PRECISION  FACT, RSUM
      INTEGER           ICOMP, ID, IR, J, JCOL, JD, L, LL, MJ
C     .. Local Arrays ..
      DOUBLE PRECISION  BASM(5), HB(7,4)
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
C Compute local basis
C
      FACT = 1.D0
      BASM(1) = 1.D0
      DO 40 L = 1, MAXORD
         FACT = FACT*H/DBLE(L)
         BASM(L+1) = FACT
         DO 20 J = 1, KCOL
            HB(J,L) = FACT*B(J,L)
   20    CONTINUE
   40 CONTINUE
C
C Ininitialize left hand half to 0 and right hand half to I
C
      DO 80 J = 1, MSTAR
         DO 60 IR = 1, MSTAR
            GI(IR,J) = 0.D0
            GI(IR,MSTAR+J) = 0.D0
   60    CONTINUE
         GI(J,MSTAR+J) = 1.D0
   80 CONTINUE
C
C Compute the block gi
C
      IR = 1
      DO 180 ICOMP = 1, NEQ
         MJ = M(ICOMP)
         IR = IR + MJ
         DO 160 L = 1, MJ
            ID = IR - L
            DO 120 JCOL = 1, MSTAR
               RSUM = 0.D0
               DO 100 J = 1, KCOL
                  RSUM = RSUM - HB(J,L)*VI(ICOMP,J,JCOL)
  100          CONTINUE
               GI(ID,JCOL) = RSUM
  120       CONTINUE
            JD = ID - 1
            DO 140 LL = 1, L
               GI(ID,JD+LL) = GI(ID,JD+LL) - BASM(LL)
  140       CONTINUE
  160    CONTINUE
  180 CONTINUE
      RETURN
      END
