      SUBROUTINE D03PFG(PLFBND,BNDPFF,T,G,U,UB,X,NPDE,NPTS,IBND,NV,V,VD,
     *                  IRES,IIFLAG)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-------------------------------------------------------------------
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  BC description routine to call external routines used in D03PFF,
C  D03PLF, D03PSF.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IBND, IIFLAG, IRES, NPDE, NPTS, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  G(NPDE), U(NPDE,NPTS), UB(NPDE,3), V(*), VD(*),
     *                  X(NPTS)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDPFF, PLFBND
C     .. Local Scalars ..
      INTEGER           J, K
C     .. Executable Statements ..
C     If using D03PFF then pass limited number of solution values to use
      IF (IBND.EQ.0) THEN
         DO 40 J = 1, 3
            DO 20 K = 1, NPDE
               UB(K,J) = U(K,J)
   20       CONTINUE
   40    CONTINUE
      ELSE
         DO 80 J = 1, 3
            DO 60 K = 1, NPDE
               UB(K,J) = U(K,NPTS-J+1)
   60       CONTINUE
   80    CONTINUE
      END IF
C
      CALL PLFBND(NPDE,NPTS,T,X,UB,IBND,G,IRES)
      CALL BNDPFF(NPDE,NPTS,T,X,U,NV,V,VD,IBND,G,IRES)
C
      IF ((IRES.LT.-1) .OR. (IRES.EQ.0) .OR. (IRES.GT.3)) THEN
         IIFLAG = 2
      END IF
C
      RETURN
      END
