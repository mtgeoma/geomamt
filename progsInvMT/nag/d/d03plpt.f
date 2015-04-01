      SUBROUTINE D03PLP(NPDE,T,X,U,UX,NV,V,VDOT,P,C,D,S,IRES)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C----------------------------------------------------------------------
C  PDEDEF routine for simple hyperbolic system, for possible use in D03P
C  and D03PSF.
C----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, X
      INTEGER           IRES, NPDE, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NPDE), D(NPDE), P(NPDE,NPDE), S(NPDE),
     *                  U(NPDE), UX(NPDE), V(*), VDOT(*)
C     .. Local Scalars ..
      INTEGER           I, J
C     .. Executable Statements ..
      DO 40 J = 1, NPDE
         C(J) = 0.0D0
         D(J) = 0.0D0
         S(J) = 0.0D0
         DO 20 I = 1, NPDE
            IF (I.EQ.J) THEN
               P(I,J) = 1.0D0
            ELSE
               P(I,J) = 0.0D0
            END IF
   20    CONTINUE
   40 CONTINUE
      RETURN
      END
