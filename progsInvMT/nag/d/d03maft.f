      SUBROUTINE D03MAF(H,M,N,NB,NPTS,PLACES,INDEX,IDIM,IN,DIST,LD,
     *                  IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     **************************************************************
C
C     H IS THE LENGTH OF THE SIDES OF THE TRIANGLES IN THE ORIGINAL
C     GRID.
C     ON THE LINE Y=-H/2 THERE ARE M GRID POINTS AT X=0,H*SQRT(3.0),
C     . . .,(M-1)*H*SQRT(3.0) IN THE ORIGINAL GRID.
C     ON THE LINE X=0 THERE ARE N GRID POINTS AT Y=-H/2,H/2,...,
C     (N-3/2)*H IN THE ORIGINAL GRID.
C     WE WILL REFER TO A LINE OF POINTS OF THE ORIGINAL GRID WITH
C     THE SAME COORDINATES AS A GRID COLUMN.
C     THE REGION OVER WHICH THE EQUATION IS TO BE SOLVED MUST LIE
C     WITHIN THE RECTANGLE WITH CORNERS AT (0,0),(0,(N-1)*H),
C     ((M-1)*H*SQRT(3.0),0).
C     NB      IS THE NUMBER OF TIMES A TRIANGLE SIDE IS BISECTED
C             TO FIND ITS POINT OF INTERSECTION WITH THE BOUNDARY.
C             THE ACCURACY IS THEREFORE H/2**NB.
C     NPTS    IS THE NUMBER OF POINTS AT WHICH THE SOLUTION IS
C             FOUND. THIS IS CALCULATED BY D03MAF. IN THE EVENT
C             OF FAILURE BECAUSE NDIM IS TOO SMALL (IFAIL=1) OR
C             BECAUSE AN INTERNAL POINT IS FOUND ON OR OUTSIDE THE
C             BOUNDING RECTANGLE (IFAIL=2) OR BECAUSE LD IS TOO
C             SMALL (IFAIL=3) THEN NPTS IS SET TO ZERO.
C     PLACES  IS AN ARRAY IN WHICH THE X AND Y COORDINATES OF THE
C             GRID POINTS ARE STORED. ITS SECOND DIMENSION SHOULD
C             BE AT LEAST IDIM WHICH MUST BE LESS THAN NPTS.
C             D03MAF CHECKS THAT IDIM.GE.NPTS.
C     FOR THE JTH POINT AND ITS THREE NEIGHBOURS AHEAD OF IT IN
C     THE ORDERING, INDEX(I,J),I=1,4 HOLDS
C             (THE POINT NUMBER) FOR INTERNAL POINTS
C             - (THE POINT NUMBER) FOR BOUNDARY POINTS AND
C             ZERO FOR EXTERNAL POINTS.
C     IN      IS THE NAME OF AN INTEGER VALUED FUNCTION TO BE
C             SUPPLIED BY THE USER. IT HAS REAL ARGUMENTS X,Y AND
C             SHOULD RETURN THE VALUE 1 IF (X,Y) LIES INSIDE THE
C             REGION AND 0 IF OUTSIDE.
C     DIST    IS AN ARRAY USED FOR WORKSPACE. ITS SECOND DIMENSION
C             SHOULD BE AT LEAST LD WHICH MUST BE AT LEAST 4*N.
C
C     **************************************************************
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D03MAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H
      INTEGER           IDIM, IFAIL, LD, M, N, NB, NPTS
C     .. Array Arguments ..
      DOUBLE PRECISION  DIST(4,LD), PLACES(2,IDIM)
      INTEGER           INDEX(4,IDIM)
C     .. Function Arguments ..
      INTEGER           IN
      EXTERNAL          IN
C     .. Local Scalars ..
      DOUBLE PRECISION  CASE, DISTM, DMIN, H2, HH, HX, RT3, TEMP, X,
     *                  XBDY, XIN, XOUT, Y, YBDY, YIN, YOUT
      INTEGER           I, I1, I2, I3, IERROR, IND0, ISAVE, IX, IXX, IY,
     *                  IYX, J, J1, K, K1, L, MIND0, MM1, MMM, N2, N3,
     *                  N4, NP
C     .. Local Arrays ..
      DOUBLE PRECISION  D(6), XSTEP(6), YSTEP(6)
      INTEGER           IND(6)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN, MOD, DBLE, SQRT
C     .. Executable Statements ..
      ISAVE = IFAIL
      NPTS = 0
C     CHECK THE PARAMETERS M,N,NB
      IERROR = 4
      IF (M.LE.2) GO TO 920
      IERROR = 5
      IF (N.LE.2) GO TO 920
      IERROR = 6
      IF (NB.LE.0) GO TO 920
C
C     SET UP SOME USEFUL PARAMETERS.
C
      RT3 = SQRT(3.0D0)
      HH = H + H
      N2 = N*2
      N3 = N*3
      N4 = N*4
      MM1 = M + M + 1
      MMM = M + M - 1
      H2 = H/2.0D0
      HX = H2*RT3
      IF (LD.LT.N4) GO TO 860
      XSTEP(1) = HX
      YSTEP(1) = -H2
      XSTEP(2) = HX
      YSTEP(2) = H2
      XSTEP(3) = 0.0D0
      YSTEP(3) = H
      XSTEP(4) = -HX
      YSTEP(4) = H2
      XSTEP(5) = -HX
      YSTEP(5) = -H2
      XSTEP(6) = 0.0D0
      YSTEP(6) = -H
C
C     WE SWEEP THE POINTS OF THE ORIGINAL GRID BY COLUMNS (I.E. IN
C     ORDER OF INCREASING X AND Y WITH THE Y COORDINATE VARYING
C     THE MOST RAPIDLY). WE KEEP INFORMATION ABOUT FOUR COLUMNS
C     OF POINTS OF THE ORIGINAL GRID. DIST(I,J), I=1,3, J=1,3*N
C     HOLDS THE DISTANCE FROM THE JTH OF THESE POINTS TO THE
C     BOUNDARY IN THE FORWARD DIRECTION GIVEN BY XSTEP(I),YSTEP(I)
C     IF THIS DOES NOT EXCEED H AND OTHERWISE THE DUMMY VALUE HH.
C     DIST(4,J), J=1,3*N HOLDS FOR THE SAME POINTS THE POINT
C     NUMBER FOR INTERNAL POINTS, - (THE POINT NUMBER) FOR
C     BOUNDARY POINTS AND ZERO FOR EXTERNAL POINTS.
C     DIST(4,J), J=3*N+1,4*N HOLDS THE VALUES OF IN FOR
C     THE POINTS OF THE FOURTH COLUMN.
C
      DO 40 J = 1, N4
         DO 20 I = 1, 4
            DIST(I,J) = 0.0D0
   20    CONTINUE
   40 CONTINUE
      DO 80 I = N2, N3
         DO 60 K = 1, 3
            DIST(K,I) = HH
   60    CONTINUE
   80 CONTINUE
      DO 840 IX = 1, MM1
         IXX = MOD(IX-1,2)
C
C        FIND WHETHER THE ELEMENTS OF THE FOURTH COLUMN ARE INSIDE
C        OR OUTSIDE THE REGION.
C
         X = HX*DBLE(IX-1)
         DO 100 IY = 1, N
            I = IN(X,H2*DBLE(IY+IY+IXX-3))
            IF (I.EQ.1 .AND. (IX.EQ.1 .OR. IX.GE.MMM .OR. IY.EQ.1 .OR.
     *          IY+IXX.GT.N)) GO TO 880
            I1 = N3 + IY
            DIST(4,I1) = I
  100    CONTINUE
         IF (IX.EQ.1) GO TO 780
C
C        FIND INFORMATION ABOUT POINTS IN THE THIRD COLUMN.
C
         X = HX*DBLE(IX-2)
         DO 340 IY = 1, N
            IYX = IY - IXX
            Y = H2*DBLE(IY+IYX-2)
C
C           IND0 HOLDS IN(X,Y) AND IND(K), K=1,3 HOLD
C           IN(X + XSTEP(K),Y + YSTEP(K)).
C
            I1 = N2 + IY
            IND0 = DIST(4,I1)
            DO 120 K = 1, 3
               D(K+3) = HH
               IND(K) = 0
  120       CONTINUE
            IF (IY.NE.N) IND(3) = DIST(4,I1+1)
            I2 = N3 + IYX
            IF (IYX.GT.0) IND(1) = DIST(4,I2)
            IF (IYX.LT.N) IND(2) = DIST(4,I2+1)
            DISTM = HH
            DO 220 K = 1, 3
               TEMP = HH
               IF (IND0.EQ.IND(K)) GO TO 200
C
C              FIND THE BOUNDARY POINT ON THE LINE JOINING (X,Y)
C              TO (X + XSTEP(K),Y + YSTEP(K)).
C
               XIN = X
               YIN = Y
               XOUT = X + XSTEP(K)
               YOUT = Y + YSTEP(K)
               IF (IND0.EQ.1) GO TO 140
               XIN = XOUT
               YIN = YOUT
               XOUT = X
               YOUT = Y
  140          DO 180 I = 1, NB
                  XBDY = (XIN+XOUT)*0.5D0
                  YBDY = (YIN+YOUT)*0.5D0
                  IF (IN(XBDY,YBDY).EQ.1) GO TO 160
                  XOUT = XBDY
                  YOUT = YBDY
                  GO TO 180
  160             XIN = XBDY
                  YIN = YBDY
  180          CONTINUE
               XBDY = (XIN+XOUT)*0.5D0
               YBDY = (YIN+YOUT)*0.5D0
               TEMP = SQRT((X-XBDY)**2+(Y-YBDY)**2)
               DISTM = MIN(DISTM,TEMP)
  200          DIST(K,I1) = TEMP
  220       CONTINUE
            TEMP = DISTM
            IF (IY.NE.1) D(6) = DIST(3,I1-1)
            I3 = N + IYX
            IF (IY.GT.IXX) D(5) = DIST(2,I3)
            IF (IYX.LT.N) D(4) = DIST(1,I3+1)
            DO 240 K = 4, 6
               IF (D(K).LE.H) DISTM = MIN(DISTM,H-D(K))
  240       CONTINUE
            CASE = 0.0D0
            IF (DISTM-H2) 260, 280, 300
  260       NPTS = NPTS + 1
            IF (NPTS.GT.IDIM) GO TO 900
            CASE = -NPTS
            GO TO 320
C
C           IF THE DISTANCE TO THE BOUNDARY IS H/2 THEN WE TAKE THE
C           DISTANCE IN THE FORWARD DIRECTION TO BE SLIGHTLY LESS.
C
  280       IF (TEMP.EQ.H2) GO TO 260
  300       IF (IND0.EQ.0) GO TO 320
            NPTS = NPTS + 1
            IF (NPTS.GT.IDIM) GO TO 900
            CASE = NPTS
  320       DIST(4,I1) = CASE
  340    CONTINUE
C
C        CHECK THAT THE LAST POINT IS NEAR THE BOUNDARY.
C
         IF (CASE.GT.0.0D0) GO TO 880
         IF (IX.EQ.2) GO TO 780
C
C        NOW FIND PLACES AND INDEX FOR POINTS IN COLUMN TWO.
C
         X = HX*DBLE(IX-3)
         DO 760 IY = 1, N
            I1 = N + IY
            IND0 = DIST(4,I1)
            IF (IND0.EQ.0) GO TO 760
            IYX = IY + IXX - 1
            Y = H2*DBLE(IY+IYX-2)
C
C           FOR K=1,2,...,6 D(K) HOLDS THE DISTANCE FROM (X,Y) TO THE
C           BOUNDARY IN THE KTH GRID DIRECTION OR HH IF THIS DISTANCE
C           EXCEEDS H. IF THE KTH NEIGHBOUR OF (X,Y) IS INTERNAL/
C           BOUNDARY/EXTERNAL THEN IND(K) HOLDS THE POINT NUMBER/
C           -(THE POINT NUMBER)/ ZERO. IF LATER WE WISH TO RESTRICT THE
C           CHOICE OF DIRECTION FOR THE MOVEMENT OF (X,Y) THEN WE
C           REPLACE ONE OR MORE OF THE D(K) BY HH.
C
            DO 360 K = 1, 6
               IND(K) = 0
               D(K) = HH
  360       CONTINUE
            I2 = N2 + IYX
            IF (IYX.EQ.0) GO TO 380
            IND(5) = DIST(4,IYX)
            D(5) = DIST(2,IYX)
            IND(1) = DIST(4,I2)
            D(1) = DIST(1,I1)
  380       IF (IYX.EQ.N) GO TO 400
            IND(4) = DIST(4,IYX+1)
            D(4) = DIST(1,IYX+1)
            IND(2) = DIST(4,I2+1)
            D(2) = DIST(2,I1)
  400       IF (IY.EQ.1) GO TO 420
            IND(6) = DIST(4,I1-1)
            D(6) = DIST(3,I1-1)
  420       IF (IY.NE.N) IND(3) = DIST(4,I1+1)
            D(3) = DIST(3,I1)
            DO 440 K = 4, 6
               IF (D(K).LE.H) D(K) = H - D(K)
               IF (D(K).EQ.H2) D(K) = H
  440       CONTINUE
            IF (IND0.LT.0) GO TO 460
            PLACES(1,IND0) = X
            PLACES(2,IND0) = Y
            GO TO 720
C
C           SEE IF THERE IS ANY CHOICE FOR THE DIRECTION OF MOVEMENT OF
C           THE CURRENT POINT TO THE BOUNDARY AND COUNT THE NUMBER
C           OF INTERNAL NEIGHBOURS.
C
  460       I = 0
            J = 0
            DO 480 K = 1, 6
               IF (D(K).LE.H2) I = I + 1
               IF (IND(K).GT.0) J = J + 1
  480       CONTINUE
C
C           OMIT ANY POINT WITHOUT ANY NEIGHBOURS INSIDE THE REGION.
C
            IF (J.GT.0) GO TO 580
            NPTS = NPTS - 1
            DIST(4,I1) = 0
            DO 500 J = 1, 3
               I3 = ABS(IND(J+3))
               IF (IND(J+3).NE.0) INDEX(J+1,I3) = 0
  500       CONTINUE
            J1 = N + IY + 1
            DO 520 J = J1, N3
               TEMP = DIST(4,J)
               IF (TEMP.GT.0.0D0) DIST(4,J) = TEMP - 1.0D0
               IF (TEMP.LT.0.0D0) DIST(4,J) = TEMP + 1.0D0
  520       CONTINUE
            J1 = N + IY - 1
            DO 560 J = 1, J1
               NP = ABS(DIST(4,J))
               IF (NP.EQ.0) GO TO 560
               DO 540 I = 2, 3
                  IF (INDEX(I,NP).GT.-IND0) INDEX(I,NP) = INDEX(I,NP) -
     *                1
                  IF (INDEX(I,NP).LT.IND0) INDEX(I,NP) = INDEX(I,NP) + 1
  540          CONTINUE
  560       CONTINUE
            GO TO 760
  580       IF (I.EQ.1) GO TO 680
C
C           AVOID MOVEMENT TOWARDS ANOTHER BOUNDARY POINT.
C
            DO 600 K = 1, 6
               IF (IND(K).GE.0 .OR. D(K).GT.H2) GO TO 600
               D(K) = HH
               I = I - 1
               IF (I.EQ.1) GO TO 680
  600       CONTINUE
C
C           AVOID MOVEMENT WHICH INCREASES THE INTERNAL ANGLE IN A
C           FORWARD POINTING TRIANGLE.
C
            K1 = 6
            DO 660 K = 1, 6
               IF (IND(K1).LE.0 .OR. IND(K).GE.0) GO TO 640
C
C              THE DIRECTIONS K-2 AND K-3 ARE TO BE AVOIDED.
C
               DO 620 J = 2, 3
                  L = K - J
                  IF (L.LE.0) L = L + 6
                  IF (D(L).GT.H2) GO TO 620
                  D(L) = HH
                  I = I - 1
                  IF (I.EQ.1) GO TO 680
  620          CONTINUE
  640          K1 = K
  660       CONTINUE
C
C           CHOOSE THE DIRECTION INVOLVING LEAST MOVEMENT.
C
  680       DMIN = HH
            DO 700 K = 1, 6
               IF (D(K).GE.DMIN) GO TO 700
               L = K
               DMIN = D(K)
  700       CONTINUE
            MIND0 = -IND0
            PLACES(1,MIND0) = X + XSTEP(L)*(DMIN/H)
            PLACES(2,MIND0) = Y + YSTEP(L)*(DMIN/H)
  720       I = ABS(IND0)
            INDEX(1,I) = IND0
            DO 740 K = 1, 3
               INDEX(K+1,I) = IND(K)
  740       CONTINUE
  760    CONTINUE
C
C        SHIFT THE NUMBERS STORED IN DIST IN PREPARATION FOR THE
C        NEXT SET OF COLUMNS.
C
  780    DO 820 I = 1, 4
            DO 800 J = 1, N3
               I1 = J + N
               DIST(I,J) = DIST(I,I1)
  800       CONTINUE
  820    CONTINUE
  840 CONTINUE
      IFAIL = 0
      GO TO 940
  860 IERROR = 3
      GO TO 920
  880 IERROR = 2
      GO TO 920
  900 IERROR = 1
  920 IFAIL = P01ABF(ISAVE,IERROR,SRNAME,0,P01REC)
      NPTS = 0
  940 RETURN
      END
