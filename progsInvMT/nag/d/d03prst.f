      SUBROUTINE D03PRS(PHFPDE,PDEPCF,X,T,U,NPDE,NPTS,UMEAN,UX,G,CL,
     *                  RFUN,V,NV,VT,SPDEF1)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C***********************************************************************
C     Routine to provide estimates of the flux for use in space remesh
C     SKFLUX routine from SPRINT
C  ---------------------------------------------------------------------
C     In order to make available to the user an estimate of the flux
C     at the end-pts for dirichlet b.c.s it uses 3-pt end-on formulae
C     material discontinuities must not be present within these points.
C
C     X:   DIMENSION X(NPTS);
C          X(1),...,X(NPTS) is a partition of (X(1),X(NPTS))
C
C     T:   the time variable;
C
C     U:   DIMENSION U(NPDE,NPTS);
C          U(I,J) is an approximation of U(I,X(J),T),I = 1,...,NPDE;
C                                                    J = 1,...,NPTS;
C     NPDE: the number of partial differential equations;
C
C     NPTS: the number of gridpoints;
C
C     All other arrays are workspace except for
C     RFUN(NPDE,NPTS+1): empty on entry, contains flux values on exit.
C            RFUN(J,I) contains the flux for the Jth PDE at
C                      I = 1 or NPTS+1 the left or the right boundary
C                      I = 2,NPTS midway between X(I) and X(I-1).
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           NPDE, NPTS, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  CL(NPDE,NPDE), G(NPDE), RFUN(NPDE,NPTS+1),
     *                  U(NPDE,NPTS), UMEAN(NPDE), UX(NPDE), V(*),
     *                  VT(*), X(NPTS)
C     .. Subroutine Arguments ..
      EXTERNAL          PDEPCF, PHFPDE, SPDEF1
C     .. Local Scalars ..
      DOUBLE PRECISION  H, HQ, HSUM, XMEAN
      INTEGER           I, J, K, KK, L
C     .. Executable Statements ..
      K = 1
      DO 60 L = 2, NPTS
         H = X(L) - X(L-1)
C        Assume non-polar case for flux calculation.
         XMEAN = (X(L)+X(L-1))*0.5D0
         IF (L.EQ.2) THEN
C           First segment
            HQ = X(3) - X(2)
            HSUM = H + HQ
            DO 20 J = 1, NPDE
C              Est. flux at left end point by 3-pt end on
C              formula  used via SPDEF1 to evaluate the flux rfun
               UX(J) = -(2.0D0*H+HQ)*U(J,1)/(H*HSUM) + HSUM*U(J,2)
     *                 /(H*HQ) - H*U(J,3)/(HQ*HSUM)
   20       CONTINUE
            CALL SPDEF1(PHFPDE,PDEPCF,T,X(1),NPDE,U,UX,CL,G,RFUN,NV,V,
     *                  VT,K,KK)
         END IF
         DO 40 J = 1, NPDE
            UMEAN(J) = (U(J,L-1)+U(J,L))*0.5D0
            UX(J) = (U(J,L)-U(J,L-1))/H
   40    CONTINUE
         CALL SPDEF1(PHFPDE,PDEPCF,T,XMEAN,NPDE,UMEAN,UX,CL,G,RFUN(1,L),
     *               NV,V,VT,K,KK)
   60 CONTINUE
C
C     Calculate flux at rightmost boundary
C
      H = X(NPTS-1) - X(NPTS-2)
      HQ = X(NPTS) - X(NPTS-1)
      HSUM = H + HQ
      DO 80 J = 1, NPDE
C        Calculate UX and the flux at the end-pt for remesh using
C        3-pt end-on formula is used via SPEDF1 to form flux RFUN
         UX(J) = (2.0D0*HQ+H)*U(J,NPTS)/(HQ*HSUM) - HSUM*U(J,NPTS-1)
     *           /(H*HQ) + HQ*U(J,NPTS-2)/(H*HSUM)
   80 CONTINUE
      I = NPTS + 1
      CALL SPDEF1(PHFPDE,PDEPCF,T,X(NPTS),NPDE,U(1,NPTS),UX,CL,G,
     *            RFUN(1,I),NV,V,VT,K,KK)
      RETURN
      END
