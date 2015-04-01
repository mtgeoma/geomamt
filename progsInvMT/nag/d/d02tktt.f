      SUBROUTINE D02TKT(X,ZVAL,A,COEF,XLEFT,XRIGHT,Z,DMZ,KCOL,NEQ,
     *                  MAXORD,M,MSTAR,GETLB,GETDM,DMVAL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C**********************************************************************
C
C   Purpose:
C                                 (1)         (m(1)-1)
C      Evaluate z(u(x))=( u (x), u (x), ... ,u (x),       .....
C                          1      1           1
C                                        (m(neq)-1)
C                               ......, u   (x)         )
C                                        neq
C      at one point x.
C
C   Arguments:
C     X      - point at which z(u(x)) required
C     ZVAL   - the computed values z(u(x)
C     A      - array of mesh independent rk-basis coefficients
C     COEF   - array of mesh dependent monomial coefficients
C     XLEFT  - left hand point of interval containing X
C     XRIGHT - right hand point of interval containing X
C     Z      - the current solution vector for the interval
C     DMZ    - the array of m(j)-th derivatives of the current
C              solution at each collocation point in the interval
C     KCOL   - the number of collocation points
C     NEQ    - the number of ODEs
C     MAXORD - the maximal order of the ODEs
C     M      - the orders of the ODEs
C     MSTAR  - the total number of unkowns (sum(m(i)),i=1,n)
C     GETLB  - determines if local rk-basis is required
C     GETDM  - determines if DMVAL should be computed
C     DMVAL  - the m(j)-th derivatives of each ODE computed at X
C
C   Author:
C      R.W. Brankin, NAG Ltd., August 1994
C      (modified version of COLNEW routine APPROX)
C
C**********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X, XLEFT, XRIGHT
      INTEGER           KCOL, MAXORD, MSTAR, NEQ
      LOGICAL           GETDM, GETLB
C     .. Array Arguments ..
      DOUBLE PRECISION  A(7,4), COEF(KCOL,KCOL), DMVAL(NEQ),
     *                  DMZ(NEQ,KCOL), Z(MSTAR), ZVAL(MSTAR)
      INTEGER           M(NEQ)
C     .. Local Scalars ..
      DOUBLE PRECISION  FACT, S, ZSUM
      INTEGER           IZ, J, JCOMP, L, LB, LL, MJ
C     .. Local Arrays ..
      DOUBLE PRECISION  BM(4), DM(7)
C     .. External Subroutines ..
      EXTERNAL          D02TKY
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
C Compute mesh independent rk-basis if required
C
      IF (GETLB) THEN
         S = (X-XLEFT)/(XRIGHT-XLEFT)
         CALL D02TKY(S,COEF,KCOL,MAXORD,A,GETDM,DM)
      END IF
C
C Compute mesh dependent rk-basis.
C
      BM(1) = X - XLEFT
      DO 20 L = 2, MAXORD
         BM(L) = BM(1)/DBLE(L)
   20 CONTINUE
C
C Evaluate  z( u(x) ).
C
      IZ = 1
      DO 100 JCOMP = 1, NEQ
         MJ = M(JCOMP)
         IZ = IZ + MJ
         DO 80 L = 1, MJ
            ZSUM = 0.D0
            DO 40 J = 1, KCOL
               ZSUM = ZSUM + A(J,L)*DMZ(JCOMP,J)
   40       CONTINUE
            DO 60 LL = 1, L
               LB = L + 1 - LL
               ZSUM = ZSUM*BM(LB) + Z(IZ-LL)
   60       CONTINUE
            ZVAL(IZ-L) = ZSUM
   80    CONTINUE
  100 CONTINUE
C
C Evaluate  dmval(j) (= m(j)-th derivative of uj) ?
C
      IF (GETDM) THEN
         DO 120 JCOMP = 1, NEQ
            DMVAL(JCOMP) = 0.D0
  120    CONTINUE
         DO 160 J = 1, KCOL
            FACT = DM(J)
            DO 140 JCOMP = 1, NEQ
               DMVAL(JCOMP) = DMVAL(JCOMP) + FACT*DMZ(JCOMP,J)
  140       CONTINUE
  160    CONTINUE
      END IF
      RETURN
      END
