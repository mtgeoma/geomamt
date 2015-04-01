      SUBROUTINE D02TKJ(N,MSTAR,KCOL,Z,XI,SCALE,DSCALE,M,NEQ,MAXORD)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C**********************************************************************
C
C   Purpose:
C      Provide a proper scaling of the state variables, used
C      to control the damping factor for a newton iteration.
C
C   arguments:
C      N      - number of mesh subintervals
C      MSTAR  - number of unknomns in z(u(x))  (sum(M(i),1=1,NEQ)
C      KCOL   - number of collocation points
C      Z      - the global unknown vector      (covers whole mesh)
C      XI     - the current mesh
C      SCALE  - scaling vector for z
C      DSCALE - scaling vector for dmz
C      M      - orders of the odes
C      NEQ    - the number of odes
C      MAXORD - maximal order of the odes
C
C   Author:
C      R.W. Brankin, NAG Ltd., August 1994
C      (modified version of COLNEW routine SKALE)
C
C**********************************************************************
C     .. Scalar Arguments ..
      INTEGER           KCOL, MAXORD, MSTAR, N, NEQ
C     .. Array Arguments ..
      DOUBLE PRECISION  DSCALE(NEQ,KCOL,N), SCALE(MSTAR,N+1), XI(N+1),
     *                  Z(MSTAR,N+1)
      INTEGER           M(NEQ)
C     .. Local Scalars ..
      DOUBLE PRECISION  H, SCAL
      INTEGER           ICOMP, IDMZ, IZ, J, L, MJ, NP1
C     .. Local Arrays ..
      DOUBLE PRECISION  BASM(5)
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Executable Statements ..
C
      BASM(1) = 1.D0
      DO 100 J = 1, N
         IZ = 1
         H = XI(J+1) - XI(J)
         DO 20 L = 1, MAXORD
            BASM(L+1) = BASM(L)*H/DBLE(L)
   20    CONTINUE
         DO 80 ICOMP = 1, NEQ
            SCAL = (ABS(Z(IZ,J))+ABS(Z(IZ,J+1)))*.5D0 + 1.D0
            MJ = M(ICOMP)
            DO 40 L = 1, MJ
               SCALE(IZ,J) = BASM(L)/SCAL
               IZ = IZ + 1
   40       CONTINUE
            SCAL = BASM(MJ+1)/SCAL
C            DO 60 IDMZ = ICOMP, KD, NEQ
C               DSCALE(IDMZ,J) = SCAL
C   60       CONTINUE
            DO 60 IDMZ = 1, KCOL
               DSCALE(ICOMP,IDMZ,J) = SCAL
   60       CONTINUE
   80    CONTINUE
  100 CONTINUE
      NP1 = N + 1
      DO 120 IZ = 1, MSTAR
         SCALE(IZ,NP1) = SCALE(IZ,N)
  120 CONTINUE
      RETURN
      END
