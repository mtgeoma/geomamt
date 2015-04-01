      SUBROUTINE D03PHK(XP,UP,IPTS,X,U,NPTS,NPDE,IFLAG)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C---------------------------------------------------------------------
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  Routine to compute flux at the interpolation points XP(IP)
C  linear interpolation is used.
C
C  Parameter list :
C
C  XP(IPTS) ;  the spatial interpolation points , for efficiency
C               XP(I) < XP(I+1) , I = 1, IPTS-1.
C  UP(NPDE,IPTS);   array to hold the values found by linear interp
C
C  X(NPTS)  ;  is an array that contains the original mesh.
C
C  U(NPDE, NPTS+1) ; the original flux ; U(J,1) are the fluxes at
C                   X(1) and U(J,NPTS) are the fluxes at X(NPTS)
C                   while  U(J,K) is the flux
C                   at 0.5 * (X(K-1) + X(K)) , K = 2, NPTS , J = 1,NPDE.
C  IFLAG  ;  error flag ; set  to 1 if extrapolation tried else 0.
C
C       IT          HAS VALUE 1 OR 2 DEPENDING ON HOW MANY COMPONENTS
C                   OF THE ARRAY UP ARE REQUIRED.
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IFLAG, IPTS, NPDE, NPTS
C     .. Array Arguments ..
      DOUBLE PRECISION  U(NPDE,*), UP(NPDE,IPTS), X(NPTS), XP(IPTS)
C     .. Scalars in Common ..
      DOUBLE PRECISION  DUNFLO, UROUND
      INTEGER           IOVFLO
C     .. Local Scalars ..
      DOUBLE PRECISION  UI, XI1, XL, XR
      INTEGER           I, J, JI
C     .. Common blocks ..
      COMMON            /CD03PC/DUNFLO, UROUND, IOVFLO
C     .. Save statement ..
      SAVE              /CD03PC/
C     .. Executable Statements ..
      DO 100 JI = 1, IPTS
         IFLAG = 0
C
C ... Search the grid points which bound XP(JI) ...
C
         IF (XP(JI).LT.X(1)-UROUND) GO TO 40
         DO 20 I = 1, NPTS
            XI1 = X(I)
            IF (I.LT.NPTS) XI1 = (XI1+X(I+1))*0.5D0
            IF (XP(JI).LE.XI1) GO TO 60
   20    CONTINUE
C
         IF (XP(JI).GT.X(NPTS)+UROUND) GO TO 40
         I = NPTS
         GO TO 60
C       XP(JI) IS OUTSIDE THE MESH
   40    CONTINUE
         IFLAG = 1
         RETURN
   60    CONTINUE
C
C ... The point XP(JI) Is bounded by the i-1th and the ith  pts ...
C
         IF (I.GT.1) THEN
            IF (I.EQ.NPTS) THEN
C
C    I=NPTS, INTERPOLATE U(J,NPTS), U(J,NPTS+1) IN (LAST MID-PT, X(NPTS)
C
               XL = 0.5D0*(X(NPTS-1)+X(NPTS))
               XR = X(NPTS)
            ELSE
C
C  I=2..NPTS, INTERPOLATE U(J,I), U(J,I+1) IN (LAST MID-PT, NEXT MID-PT)
C
               XL = 0.5D0*(X(I)+X(I-1))
               XR = XI1
            END IF
         ELSE
C
C          I=1, INTERPOLATE U(J,1) AND U(J,2) IN (X(1),0.5(X(1)+X(2)))
C
            XL = X(1)
            XR = XI1
         END IF
         DO 80 J = 1, NPDE
            UI = U(J,I)
C             APPLY LINEAR INTERPOLATION
            UP(J,JI) = UI + (XP(JI)-XL)*(U(J,I+1)-UI)/(XR-XL)
   80    CONTINUE
  100 CONTINUE
  120 CONTINUE
      RETURN
      END
