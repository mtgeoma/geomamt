      SUBROUTINE D03PCW(X,NPTS,NC,APLUS,BPLUS,EPLUS,XHAT,PR,DPL,DPR)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C----------------------------------------------------------------------
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   Routine to define the coefficients associated with the spatial mesh
C   of the October 1986 NAG Library  D03 Chapter
C   implementation of the new Skeel (1986) discretisation.
C   This routine must be recalled whenever the spatial mesh is changed.
C
C     INPUT:
C
C     X    ; DIMENSION X(NPTS);
C            X(1),...,X(NPTS) is a partition of (X(1),X(NPTS))
C
C     NPTS ; the number of gridpoints
C
C     NC   ;  a number designing the kind of space coordinates;
C             NC = 0 : Cartesian coordinates;
C             NC = 1 : polar coordinates;
C             NC = 2 : spherical coordinates;
C
C     Output the coeffs associated with the spatial mesh - the arrays:
C
C     APLUS(NPTS), BPLUS(NPTS), EPLUS(NPTS), XHAT(NPTS), PR(NPTS)
C     DPL(NPTS), DPR(NPTS)
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           NC, NPTS
C     .. Array Arguments ..
      DOUBLE PRECISION  APLUS(NPTS), BPLUS(NPTS), DPL(NPTS), DPR(NPTS),
     *                  EPLUS(NPTS), PR(NPTS), X(NPTS), XHAT(NPTS)
C     .. Local Scalars ..
      DOUBLE PRECISION  H, XL, XR, XSUM
      INTEGER           L
      LOGICAL           SINGLR
C     .. Intrinsic Functions ..
      INTRINSIC         LOG
C     .. Executable Statements ..
C
      SINGLR = .FALSE.
      IF (NC.GT.0 .AND. X(1).LE.1.0D-8) SINGLR = .TRUE.
C
      APLUS(1) = 0.0D0
      BPLUS(1) = 0.0D0
      EPLUS(1) = 0.0D0
      PR(1) = 0.0D0
      DPR(1) = 0.0D0
      DPL(1) = 0.0D0
      XHAT(1) = 0.0D0
C
      IF (NC.EQ.0) THEN
C
C ... Non-polar case ...
C
         DO 20 L = 2, NPTS
            H = X(L) - X(L-1)
            XSUM = X(L) + X(L-1)
            APLUS(L) = H*0.5D0
            BPLUS(L) = H*0.5D0
            EPLUS(L) = 1.0D0
            PR(L) = 0.5D0
            DPR(L) = 1.D0/H
            DPL(L) = -DPR(L)
            XHAT(L) = XSUM*0.5D0
   20    CONTINUE
C
         RETURN
      END IF
C
C ... Polar singularity assume that X(1) = 0.0D0 ...
C
      IF (SINGLR) THEN
         XHAT(2) = X(2)*2.D0/3.D0
         EPLUS(2) = 1.5D0/X(2)
         APLUS(2) = 1.0D0/(NC+1)
         BPLUS(2) = X(2)**(NC+1)/(NC+1)
         PR(2) = 4.0D0/9.0D0
         DPR(2) = 2*XHAT(2)/(X(2)**2)
         DPL(2) = -DPR(2)
      END IF
C
C ... Cylindrical polar co-ordinates ...
C
      IF (NC.EQ.1) THEN
         IF (SINGLR) THEN
C
C ... M = 1  case polar singularity present ...
C
            DO 40 L = 3, NPTS
               XR = X(L)
               XL = X(L-1)
               XSUM = X(L) + X(L-1)
               H = X(L) - X(L-1)
               DPR(L) = 1.D0/LOG(X(L)/X(L-1))
               APLUS(L) = -XL**2*0.5D0 + XSUM*H*0.25D0*DPR(L)
               BPLUS(L) = XR**2*0.5D0 - XSUM*H*0.25D0*DPR(L)
               XHAT(L) = 2.0D0*(XR**2+XR*XL+XL**2)/(3.0D0*XSUM)
               EPLUS(L) = XSUM*H*0.5D0*DPR(L)/XHAT(L)
               PR(L) = (XHAT(L)**2-XL**2)/(XSUM*H)
               DPR(L) = 2*XHAT(L)/(XSUM*H)
               DPL(L) = -DPR(L)
   40       CONTINUE
C
         ELSE
C
C ... M = 1  case no polar singularity present ...
C
            DO 60 L = 2, NPTS
               XR = X(L)
               XL = X(L-1)
               XSUM = X(L) + X(L-1)
               H = X(L) - X(L-1)
               DPR(L) = 1.D0/LOG(X(L)/X(L-1))
               APLUS(L) = -XL**2*0.5D0 + XSUM*H*0.25D0*DPR(L)
               BPLUS(L) = XR**2*0.5D0 - XSUM*H*0.25D0*DPR(L)
               XHAT(L) = H*DPR(L)
               EPLUS(L) = XHAT(L)
               PR(L) = LOG(XHAT(L)/X(L-1))*DPR(L)
               DPR(L) = 1.D0/(X(L)-X(L-1))
               DPL(L) = -DPR(L)
   60       CONTINUE
C
         END IF
C
C ... Spherical polar co-ordinates ...
C
      ELSE IF (NC.EQ.2) THEN
         IF (SINGLR) THEN
C
C ... M = 2 case singularity present ...
C
            DO 80 L = 3, NPTS
               XR = X(L)
               XL = X(L-1)
               XSUM = X(L) + X(L-1)
               H = X(L) - X(L-1)
               APLUS(L) = -XL/6.0D0*(2.0D0*XL**2-XL*XR-XR**2)
               BPLUS(L) = XR/6.0D0*(2.0D0*XR**2-XL*XR-XL**2)
               XHAT(L) = 2.0D0*(XR**2+XR*XL+XL**2)/(3.0D0*XSUM)
               EPLUS(L) = XSUM*XL*XR*0.5D0/XHAT(L)
               PR(L) = (XHAT(L)**2-XL**2)/(XSUM*H)
               DPR(L) = 2*XHAT(L)/(XSUM*H)
               DPL(L) = -DPR(L)
   80       CONTINUE
C
         ELSE
C
C ... M = 2, no singularity present case ...
C
            DO 100 L = 2, NPTS
               XR = X(L)
               XL = X(L-1)
               XSUM = X(L) + X(L-1)
               H = X(L) - X(L-1)
               APLUS(L) = -XL/6.0D0*(2.0D0*XL**2-XL*XR-XR**2)
               BPLUS(L) = XR/6.0D0*(2.0D0*XR**2-XL*XR-XL**2)
               XHAT(L) = LOG(XR/XL)/H*XR*XL
               EPLUS(L) = XHAT(L)**2
               PR(L) = XR*(1.D0-XL/XHAT(L))/H
               DPR(L) = XL/(XHAT(L)**2*H)
               DPL(L) = XR*DPR(L)/XL
  100       CONTINUE
C
         END IF
C
      END IF
C
      RETURN
      END
