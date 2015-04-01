      SUBROUTINE D03PLK(NPDE,NPTS,X,T,U,UMEAN,ULEFT,URIGHT,NV,V,FLXPFF,
     *                  FLXPLF,RFLUX,IR,IRES)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-----------------------------------------------------------------------
C
C   Routine to calculate convective flux components for the modified
C   form of the Skeel discretisation for hyperbolic PDEs. Uses Van Leer
C   interpolation on the solution and then calculates the flux at the
C   mid points by calling the user-supplied flux routine.
C
C     INPUT:
C     X:   dimension X(NPTS);
C          X(1),...,X(NPTS) is a partition of (X(1),X(NPTS))
C     T:   the time variable;
C     U:   dimension U(NPDE,NPTS);
C          U(I,J) is an approximation of U(I,X(J),T),I = 1,...,NPDE;
C                                                    J = 1,...,NPTS;
C     NPDE:   the number of partial differential equations;
C     NV:   the number of coupled odes;
C     NPTS:   the number of gridpoints;
C     FLXPFF,FLXPLF: user-supplied subroutine and dummy to calculate
C                    the convective flux (RFLUX) given ULEFT and URIGHT
C                    at each midpoint.
C     OUTPUT:
C     ULEFT:  dimension ULEFT(NPDE);
C             calculated values of U at left side of midpoint;
C     URIGHT: dimension URIGHT(NPDE);
C             calculated values of U at right side of midpoint;
C             (ULEFT and URIGHT are required by routine FLXPFF/PLF)
C     RFLUX:  dimension (NPDE,NPTS); contains the values of the flux.
C             RFLUX(J,I) contains the flux for the Jth P.D.E. at
C                      I = 2,NPTS midway between X(I) and X(I-1).
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IR, IRES, NPDE, NPTS, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  RFLUX(NPDE,NPTS), U(NPDE,NPTS), ULEFT(NPDE),
     *                  UMEAN(NPDE), URIGHT(NPDE), V(*), X(NPTS)
C     .. Subroutine Arguments ..
      EXTERNAL          FLXPFF, FLXPLF
C     .. Scalars in Common ..
      INTEGER           IDEV, IIFLAG, ITRACE
C     .. Local Scalars ..
      DOUBLE PRECISION  H, HL, PH, PHIM1, PL, PR, T1, TOL, XL, XMEAN, XR
      INTEGER           I, IZ, J, K0, K1, K2, L, L1
C     .. External Subroutines ..
      EXTERNAL          D03PFR
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AD02NM/ITRACE, IDEV
      COMMON            /XD03PC/IIFLAG
C     .. Save statement ..
      SAVE              /AD02NM/, /XD03PC/
C     .. Executable Statements ..
      IF (IRES.EQ.-1) THEN
         DO 40 I = 1, NPTS
            DO 20 J = 1, NPDE
               RFLUX(J,I) = 0.0D0
   20       CONTINUE
   40    CONTINUE
         RETURN
      END IF
      IZ = 1
      TOL = 0.1D-15
C      TOL = 0.1D-10
      XR = X(1)
      DO 80 L = 2, NPTS
         XL = XR
         L1 = L - 1
         XR = X(L)
         H = XR - XL
         PR = 0.5D0
         PL = 1.0D0 - PR
         XMEAN = XL + H*PR
C      call flux routine midway between points X(L) and X(L-1) ..
         DO 60 J = 1, NPDE
            UMEAN(J) = PL*U(J,L1) + PR*U(J,L)
C      calculate left values at midpoint ..
            IF (L.EQ.2) THEN
C               leftmost point
               PHIM1 = U(J,2) - U(J,1)
               ULEFT(J) = U(J,1)
               ULEFT(J) = ULEFT(J) + PHIM1*0.5D0
            ELSE
C               upwind differencing towards the left
               K2 = L - 2
               K1 = L1
               K0 = L
               PH = U(J,K0) - U(J,K1)
               PHIM1 = U(J,K1) - U(J,K2)
               ULEFT(J) = U(J,K1)
               IF (ABS(PH).GT.TOL .AND. ABS(PHIM1).GT.TOL) THEN
                  T1 = ABS(PH/PHIM1)
                  H = X(K0) - X(K1)
                  HL = H/(X(K1)-X(K2))
                  ULEFT(J) = ULEFT(J) + PHIM1*HL*0.5D0*(PH/PHIM1+T1)
     *                       /(HL+T1)
               END IF
            END IF
C    calculate right values at midpoint
            IF (L.EQ.NPTS) THEN
C               rightmost point
               PHIM1 = U(J,NPTS) - U(J,NPTS-1)
               URIGHT(J) = U(J,NPTS-1)
               URIGHT(J) = URIGHT(J) + PHIM1*0.5D0
            ELSE
C               upwind differencing towards the right
               K2 = L + 1
               K1 = L
               K0 = L - 1
               PH = U(J,K0) - U(J,K1)
               PHIM1 = U(J,K1) - U(J,K2)
               URIGHT(J) = U(J,L)
               IF (ABS(PH).GT.TOL .AND. ABS(PHIM1).GT.TOL) THEN
                  T1 = ABS(PH/PHIM1)
                  H = X(K0) - X(K1)
                  HL = H/(X(K1)-X(K2))
                  URIGHT(J) = URIGHT(J) + PHIM1*HL*0.5D0*(PH/PHIM1+T1)
     *                        /(HL+T1)
               END IF
            END IF
   60    CONTINUE
         CALL D03PFR(FLXPFF,FLXPLF,NPDE,T,XMEAN,NV,V,ULEFT,URIGHT,
     *               RFLUX(1,L),IZ,IIFLAG)
         IF (IIFLAG.EQ.2) THEN
            IRES = IZ
            RETURN
         END IF
         IF (IZ.NE.1) THEN
            IR = IZ
            RETURN
         END IF
C        IF(ITRACE.GE.2)THEN
C           DO 50 J = 1,NPDE
C              WRITE(IDEV,139)L,ULEFT(J), URIGHT(J), RFLUX(J,L)
C 50           FORMAT(' L=',I4,' ULEFT AND RIGHT=',2d11.3,' FLUX',
C    *                D11.3)
C           CONTINUE
C        END IF
   80 CONTINUE
      RETURN
      END
