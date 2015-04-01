      SUBROUTINE D02TKN(GI,NBC,ZVAL,DGZ,GJAC,M,MAXORD,NEQ,MSTAR,DG,
     *                  DGDOTZ)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C**********************************************************************
C
C   Purpose:
C      Construct a collocation matrix row for a group of boundary
C      conditions
C
C   Arguments:
C      GI     - the sub-block of the global bvp matrix in
C               which the equations are to be formed, ie. first or last
C      NBC    - number of rows in GI (ie. no. of boundary conditions)
C      ZVAL   - z(xi)
C      DGZ    - dg . z(xi)
C      GJAC  - procedure for evaluating Jacobian of boundary conditions
C      M      - orders of ODEs
C      MAXORD - maximal order of ODEs
C      NEQ    - number of ODEs
C      DG     - Jacobian of selected boundary conditions
C      MSTAR  - number of unknowns, sum(M(i),i=1,NEQ)
C      DGDOTZ - indicates if dg . z(xi) should be evaluated
C
C   Author:
C      R.W. Brankin, NAG Ltd., August 1994
C      (modified version of COLNEW routine GDERIV)
C
C**********************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           MAXORD, MSTAR, NBC, NEQ
      LOGICAL           DGDOTZ
C     .. Array Arguments ..
      DOUBLE PRECISION  DG(NBC,NEQ,MAXORD), DGZ(NBC), GI(NBC,MSTAR),
     *                  ZVAL(NEQ,MAXORD)
      INTEGER           M(NEQ)
C     .. Subroutine Arguments ..
      EXTERNAL          GJAC
C     .. Local Scalars ..
      DOUBLE PRECISION  DOT
      INTEGER           I, J, JCOL, K
C     .. Executable Statements ..
C
      DO 60 I = 1, MAXORD
         DO 40 J = 1, NEQ
            DO 20 K = 1, NBC
               DG(K,J,I) = 0.0D0
   20       CONTINUE
   40    CONTINUE
   60 CONTINUE
C
C Evaluate jacobian dg
C
      CALL GJAC(ZVAL,NEQ,M,NBC,DG)
C
C Evaluate  dgz = dg * zval  once for a new mesh
C
      IF (DGDOTZ) THEN
         DO 120 K = 1, NBC
            DOT = 0.0D0
            DO 100 I = 1, NEQ
               DO 80 J = 1, M(I)
                  DOT = DOT + DG(K,I,J)*ZVAL(I,J)
   80          CONTINUE
  100       CONTINUE
            DGZ(K) = DOT
  120    CONTINUE
      END IF
C
C Set appropriate row of the block
C
      DO 180 K = 1, NBC
         JCOL = 0
         DO 160 I = 1, NEQ
            DO 140 J = 1, M(I)
               JCOL = JCOL + 1
               GI(K,JCOL) = DG(K,I,J)
  140       CONTINUE
  160    CONTINUE
  180 CONTINUE
C
      RETURN
      END
