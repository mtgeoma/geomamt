      SUBROUTINE D03PJG(NPDE,T,X,NPTL,U,DUDX,C,F,R,IRES)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ----------------------------------------------------------------
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Dummy routine.
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     ----------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IRES, NPDE, NPTL
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NPDE,NPDE,NPTL), DUDX(NPDE,NPTL),
     *                  F(NPDE,NPTL), R(NPDE,NPTL), U(NPDE,NPTL),
     *                  X(NPTL)
C     .. Executable Statements ..
      RETURN
C
      END
