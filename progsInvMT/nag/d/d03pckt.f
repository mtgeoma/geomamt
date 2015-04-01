      SUBROUTINE D03PCK(NPDE,T,NV,V,VDOT,NXI,XI,U,UX,RI,UTI,UTXI,VRES,
     *                  IRES)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C----------------------------------------------------------------------
C  Dummy routine to provide residual of coupled ODE system
C----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IRES, NPDE, NV, NXI
C     .. Array Arguments ..
      DOUBLE PRECISION  RI(NPDE,*), U(NPDE,*), UTI(NPDE,*),
     *                  UTXI(NPDE,*), UX(NPDE,*), V(*), VDOT(*),
     *                  VRES(*), XI(*)
C     .. Executable Statements ..
      RETURN
      END
