      SUBROUTINE D02TKQ(H,B,WI,IPVTW,KD,RHSZ,RHSDMZ,M,MAXORD,NEQ,KCOL,
     *                  MSTAR)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C**********************************************************************
C
C   Purpose:
C
C      Construct piece of condensed rhs for one subinterval.
C
C   Arguments:
C
C      H      - local stepsize.
C      B      - rk-basis coefficients
C      WI     - the sub-block of noncondensed collocation equations,
C               left-hand side part (already factorized).
C      IPVTW  - pivot info of WI factorization
C      KD     - the dimension of WI= (NEQ*KCOL)
C      RHSDMZ - the inhomogenous term of the uncondensed collocation
C               equations (as 2-d array for clarity)
C      RHSZ   - the inhomogenous term of the condensed collocation
C               equations.
C      M      - orders of ODEs
C      MAXORD - maximal order of ODEs
C      NEQ    - number of ODEs
C      KCOL   - number of collocation points (per interval)
C      MSTAR  - sum(M(i),i=1,NEQ)
C
C   Author:
C      R.W. Brankin, NAG Ltd., August 1994
C      (modified version of part of COLNEW routine GBLOCK)
C
C**********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H
      INTEGER           KCOL, KD, MAXORD, MSTAR, NEQ
C     .. Array Arguments ..
      DOUBLE PRECISION  B(7,4), RHSDMZ(NEQ,KCOL), RHSZ(MSTAR), WI(KD,KD)
      INTEGER           IPVTW(KD), M(NEQ)
C     .. Local Scalars ..
      DOUBLE PRECISION  FACT, RSUM
      INTEGER           INFO, IR, J, JCOMP, L, MJ
C     .. Local Arrays ..
      DOUBLE PRECISION  HB(7,4)
C     .. External Subroutines ..
      EXTERNAL          DGETRS
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
C Compute local basis
C
      FACT = 1.D0
      DO 40 L = 1, MAXORD
         FACT = FACT*H/DBLE(L)
         DO 20 J = 1, KCOL
            HB(J,L) = FACT*B(J,L)
   20    CONTINUE
   40 CONTINUE
C
C Note: RHSDMZ passed as 1-d array of length KD (=NEQ*KCOL)
C
      CALL DGETRS('n',KD,1,WI,KD,IPVTW,RHSDMZ,KD,INFO)
C
      IR = 1
      DO 100 JCOMP = 1, NEQ
         MJ = M(JCOMP)
         IR = IR + MJ
         DO 80 L = 1, MJ
            RSUM = 0.D0
            DO 60 J = 1, KCOL
               RSUM = RSUM + HB(J,L)*RHSDMZ(JCOMP,J)
   60       CONTINUE
            RHSZ(IR-L) = RSUM
   80    CONTINUE
  100 CONTINUE
      RETURN
      END
