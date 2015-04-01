      SUBROUTINE D03EEZ(NGX,NGY,LDA,A,RHS,XMIN,XMAX,YMIN,YMAX,PDEF,BNDY,
     *                  SCHEME,IERROR,IFAIL1)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO, HALF
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0)
      INTEGER           BOTTOM, RIGHT, TOP, LEFT
      PARAMETER         (BOTTOM=0,RIGHT=1,TOP=2,LEFT=3)
      INTEGER           S, SE, W, O, E, NW, N
      PARAMETER         (S=1,SE=2,W=3,O=4,E=5,NW=6,N=7)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XMAX, XMIN, YMAX, YMIN
      INTEGER           IERROR, IFAIL1, LDA, NGX, NGY
      CHARACTER*1       SCHEME
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,7), RHS(LDA)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDY, PDEF
C     .. Local Scalars ..
      DOUBLE PRECISION  ABND, ABNDX, ABNDY, ALPHA, BBND, BBNDX, BBNDY,
     *                  BETA, CBND, CBNDX, CBNDY, D4, D5, DELTA, EPS,
     *                  EPSLON, FF, FUDGE, GAMMA, HX, HXINV, HY, HYINV,
     *                  PHI, SCALE, SNEG, SUM, ULEFT, ULOWER, URIGHT,
     *                  USOL, UUPPER, X, XLEN, Y, YLEN
      INTEGER           I, INE, INW, ISE, ISW, J, K, KLEFT, KLOWER,
     *                  KRIGHT, KUPPER, NERR
      LOGICAL           DIADOM, ELLPTC, NEUMAN
      CHARACTER*80      REC
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN, DBLE, SIGN
C     .. Executable Statements ..
      CALL X04AAF(0,NERR)
      EPS = X02AJF()
      ELLPTC = .TRUE.
      NEUMAN = .TRUE.
      IERROR = 0
      XLEN = XMAX - XMIN
      YLEN = YMAX - YMIN
      HX = XLEN/DBLE(NGX-1)
      HY = YLEN/DBLE(NGY-1)
      HXINV = 1.0D0/HX
      HYINV = 1.0D0/HY
C
C     Set up molecule for all interior and boundary points --
C
      SNEG = 100.0D0
      DO 40 J = 1, NGY
         Y = YMIN + DBLE(J-1)*HY
         DO 20 I = 1, NGX
            X = XMIN + DBLE(I-1)*HX
            K = (J-1)*NGX + I
            CALL PDEF(X,Y,ALPHA,BETA,GAMMA,DELTA,EPSLON,PHI,FF)
            IF (PHI.NE.ZERO) NEUMAN = .FALSE.
            IF (BETA**2.GT.4.0D0*ALPHA*GAMMA) ELLPTC = .FALSE.
            IF (SCHEME.EQ.'U' .OR. SCHEME.EQ.'u') THEN
               D4 = SIGN(ONE,DELTA)
               D5 = SIGN(ONE,EPSLON)
            ELSE IF (SCHEME.EQ.'C' .OR. SCHEME.EQ.'c') THEN
               D4 = ZERO
               D5 = ZERO
            END IF
            A(K,S) = HALF*BETA/HX/HY + GAMMA/HY**2 +
     *               HALF*EPSLON*(D5-ONE)/HY
            A(K,SE) = -HALF*BETA/HX/HY
            A(K,W) = ALPHA/HX**2 + HALF*BETA/HX/HY + HALF*DELTA*(D4-ONE)
     *               /HX
            A(K,O) = -TWO*ALPHA/HX**2 - BETA/HX/HY - TWO*GAMMA/HY**2 +
     *               PHI - DELTA*D4/HX - EPSLON*D5/HY
            IF (A(K,O).LT.SNEG) SNEG = A(K,O)
            A(K,E) = ALPHA/HX**2 + HALF*BETA/HX/HY + HALF*DELTA*(D4+ONE)
     *               /HX
            A(K,NW) = -HALF*BETA/HX/HY
            A(K,N) = HALF*BETA/HX/HY + GAMMA/HY**2 +
     *               HALF*EPSLON*(D5+ONE)/HY
            RHS(K) = FF
   20    CONTINUE
   40 CONTINUE
C
C     Use automatic scaling --
C
      SCALE = -2.0D0*(1.0D0/HX**2+1.0D0/HY**2)
      SCALE = MIN(SNEG,SCALE)
C
C     Now treat the boundary points excluding the four corners.
C
C     First the top and bottom boundaries --
C
      DO 100 I = 2, NGX - 1
         X = XMIN + DBLE(I-1)*HX
         Y = YMIN
C        -- BOTTOM boundary --
         KLOWER = I
         CALL BNDY(X,Y,ABND,BBND,CBND,BOTTOM)
         IF (ABND.EQ.ZERO .AND. BBND.EQ.ZERO) THEN
            IF (IFAIL1.NE.1) THEN
               CALL X04BAF(NERR,' ** Null Boundary Condition')
               WRITE (REC,FMT='(A,1P,2E12.5)')
     *           ' ** Bottom boundary, (x,y) = ', X, Y
               CALL X04BAF(NERR,REC)
            END IF
            IERROR = 3
            RETURN
         ELSE IF (BBND.EQ.ZERO) THEN
C           -- Dirichlet Boundary Condition
C           -- Put in a trivial but scaled equation
            DO 60 J = 1, 7
               A(KLOWER,J) = ZERO
   60       CONTINUE
            A(KLOWER,O) = SCALE
            ULOWER = CBND/ABND
            RHS(KLOWER) = ULOWER*SCALE
            NEUMAN = .FALSE.
         ELSE IF (BBND.NE.ZERO) THEN
C           -- Neumann or mixed Boundary Condition
C           -- Eliminate the South point using the BC's
            A(KLOWER,N) = A(KLOWER,N) + A(KLOWER,S)
            A(KLOWER,O) = A(KLOWER,O) - TWO*HY*ABND*A(KLOWER,S)/BBND
            RHS(KLOWER) = RHS(KLOWER) - TWO*HY*CBND*A(KLOWER,S)/BBND
            A(KLOWER,S) = ZERO
            IF (ABND.NE.ZERO) NEUMAN = .FALSE.
            IF (A(KLOWER,SE).NE.ZERO .OR. A(KLOWER,NW).NE.ZERO) THEN
               IF (IFAIL1.NE.1) THEN
                  REC = ' ** Mixed Derivative in Equation '//
     *                  'and Derivative in Boundary Condition'
                  CALL X04BAF(NERR,REC)
                  WRITE (REC,FMT='(A,1P,2E12.5)')
     *              ' ** Bottom boundary, (x,y) = ', X, Y
                  CALL X04BAF(NERR,REC)
               END IF
               IERROR = 2
               RETURN
            END IF
         END IF
C        -- TOP boundary
         Y = YMAX
         KUPPER = (NGY-1)*NGX + I
         CALL BNDY(X,Y,ABND,BBND,CBND,TOP)
         IF (ABND.EQ.ZERO .AND. BBND.EQ.ZERO) THEN
            IF (IFAIL1.NE.1) THEN
               CALL X04BAF(NERR,' ** Null Boundary Condition')
               WRITE (REC,FMT='(A,1P,2E12.5)')
     *           ' ** Top boundary, (x,y) = ', X, Y
               CALL X04BAF(NERR,REC)
            END IF
            IERROR = 3
            RETURN
         ELSE IF (BBND.EQ.ZERO) THEN
C           -- Dirichlet Boundary Condition
C           -- Put in a trivial but scaled equation
            DO 80 J = 1, 7
               A(KUPPER,J) = ZERO
   80       CONTINUE
            A(KUPPER,O) = SCALE
            UUPPER = CBND/ABND
            RHS(KUPPER) = UUPPER*SCALE
            NEUMAN = .FALSE.
         ELSE IF (BBND.NE.ZERO) THEN
C           -- Neumann or mixed Boundary Condition
C           -- Eliminate the North point using the BC's
            A(KUPPER,S) = A(KUPPER,S) + A(KUPPER,N)
            A(KUPPER,O) = A(KUPPER,O) - TWO*HY*ABND*A(KUPPER,N)/BBND
            RHS(KUPPER) = RHS(KUPPER) - TWO*HY*CBND*A(KUPPER,N)/BBND
            A(KUPPER,N) = ZERO
            IF (ABND.NE.ZERO) NEUMAN = .FALSE.
            IF (A(KUPPER,NW).NE.ZERO .OR. A(KUPPER,SE).NE.ZERO) THEN
               IF (IFAIL1.NE.1) THEN
                  REC = ' ** Mixed Derivative in Equation '//
     *                  'and Derivative in Boundary Condition'
                  CALL X04BAF(NERR,REC)
                  WRITE (REC,FMT='(A,1P,2E12.5)')
     *              ' ** Top boundary, (x,y) = ', X, Y
                  CALL X04BAF(NERR,REC)
               END IF
               IERROR = 2
               RETURN
            END IF
         END IF
  100 CONTINUE
C
C     Now the left and right boundaries --
C
      DO 160 I = 2, NGY - 1
         Y = YMIN + DBLE(I-1)*HY
C        -- LEFT boundary --
         X = XMIN
         KLEFT = (I-1)*NGX + 1
         CALL BNDY(X,Y,ABND,BBND,CBND,LEFT)
         IF (ABND.EQ.ZERO .AND. BBND.EQ.ZERO) THEN
            IF (IFAIL1.NE.1) THEN
               CALL X04BAF(NERR,' ** Null Boundary Condition')
               WRITE (REC,FMT='(A,1P,2E12.5)')
     *           ' ** Left boundary, (x,y) = ', X, Y
               CALL X04BAF(NERR,REC)
            END IF
            IERROR = 3
            RETURN
         ELSE IF (BBND.EQ.ZERO) THEN
C           -- Dirichlet Boundary Condition
C           -- Put in a trivial but scaled equation
            DO 120 J = 1, 7
               A(KLEFT,J) = ZERO
  120       CONTINUE
            A(KLEFT,O) = SCALE
            ULEFT = CBND/ABND
            RHS(KLEFT) = ULEFT*SCALE
            NEUMAN = .FALSE.
         ELSE IF (BBND.NE.ZERO) THEN
C           -- Neumann or mixed Boundary Condition
C           -- Eliminate the West point using the BC's
            A(KLEFT,E) = A(KLEFT,E) + A(KLEFT,W)
            A(KLEFT,O) = A(KLEFT,O) - TWO*HX*ABND*A(KLEFT,W)/BBND
            RHS(KLEFT) = RHS(KLEFT) - TWO*HX*CBND*A(KLEFT,W)/BBND
            A(KLEFT,W) = ZERO
            IF (ABND.NE.ZERO) NEUMAN = .FALSE.
            IF (A(KLEFT,NW).NE.ZERO .OR. A(KLEFT,SE).NE.ZERO) THEN
               IF (IFAIL1.NE.1) THEN
                  REC = ' ** Mixed Derivative in Equation '//
     *                  'and Derivative in Boundary Condition'
                  CALL X04BAF(NERR,REC)
                  WRITE (REC,FMT='(A,1P,2E12.5)')
     *              ' ** Left boundary, (x,y) = ', X, Y
                  CALL X04BAF(NERR,REC)
               END IF
               IERROR = 2
               RETURN
            END IF
         END IF
C        -- RIGHT boundary
         X = XMAX
         KRIGHT = I*NGX
         CALL BNDY(X,Y,ABND,BBND,CBND,RIGHT)
         IF (ABND.EQ.ZERO .AND. BBND.EQ.ZERO) THEN
            IF (IFAIL1.NE.1) THEN
               CALL X04BAF(NERR,' ** Null Boundary Condition')
               WRITE (REC,FMT='(A,1P,2E12.5)')
     *           ' ** Right boundary, (x,y) = ', X, Y
               CALL X04BAF(NERR,REC)
            END IF
            IERROR = 3
            RETURN
         ELSE IF (BBND.EQ.ZERO) THEN
C           -- Dirichlet Boundary Condition
C           -- Put in a trivial but scaled equation
            DO 140 J = 1, 7
               A(KRIGHT,J) = ZERO
  140       CONTINUE
            A(KRIGHT,O) = SCALE
            URIGHT = CBND/ABND
            RHS(KRIGHT) = URIGHT*SCALE
            NEUMAN = .FALSE.
         ELSE IF (BBND.NE.ZERO) THEN
C           -- Neumann or mixed Boundary Condition
C           -- Eliminate the East point using the BC's
            A(KRIGHT,W) = A(KRIGHT,W) + A(KRIGHT,E)
            A(KRIGHT,O) = A(KRIGHT,O) - TWO*HX*ABND*A(KRIGHT,E)/BBND
            RHS(KRIGHT) = RHS(KRIGHT) - TWO*HX*CBND*A(KRIGHT,E)/BBND
            A(KRIGHT,E) = ZERO
            IF (ABND.NE.ZERO) NEUMAN = .FALSE.
            IF (A(KRIGHT,SE).NE.ZERO .OR. A(KRIGHT,NW).NE.ZERO) THEN
               IF (IFAIL1.NE.1) THEN
                  REC = ' ** Mixed Derivative in Equation '//
     *                  'and Derivative in Boundary Condition'
                  CALL X04BAF(NERR,REC)
                  WRITE (REC,FMT='(A,1P,2E12.5)')
     *              ' ** Right boundary, (x,y) = ', X, Y
                  CALL X04BAF(NERR,REC)
               END IF
               IERROR = 2
               RETURN
            END IF
         END IF
  160 CONTINUE
C
C     Now for the four corners--
C
C     If one of the boundary conditions at the corner is Dirichlet,
C     then use it. Otherwise use both derivatives to eliminate the
C     unknown external points.
C
      ISW = 1
      ISE = NGX
      INW = (NGY-1)*NGX + 1
      INE = NGX*NGY
C
C     Bottom Right Corner (SE) --
C
      X = XMAX
      Y = YMIN
      CALL BNDY(X,Y,ABNDX,BBNDX,CBNDX,BOTTOM)
      IF (ABNDX.EQ.ZERO .AND. BBNDX.EQ.ZERO) THEN
         IF (IFAIL1.NE.1) THEN
            CALL X04BAF(NERR,' ** Null Boundary Condition')
            WRITE (REC,FMT='(A,1P,2E12.5)')
     *        ' ** Bottom Boundary, Right End, (x,y) = ', X, Y
            CALL X04BAF(NERR,REC)
         END IF
         IERROR = 3
         RETURN
      END IF
      CALL BNDY(X,Y,ABNDY,BBNDY,CBNDY,RIGHT)
      IF (ABNDY.EQ.ZERO .AND. BBNDY.EQ.ZERO) THEN
         IF (IFAIL1.NE.1) THEN
            CALL X04BAF(NERR,' ** Null Boundary Condition')
            WRITE (REC,FMT='(A,1P,2E12.5)')
     *        ' ** Right boundary, Bottom End, (x,y) = ', X, Y
            CALL X04BAF(NERR,REC)
         END IF
         IERROR = 3
         RETURN
      END IF
      IF (BBNDX.EQ.ZERO .OR. BBNDY.EQ.ZERO) THEN
C        -- One of the Boundary Conditions is Dirichlet
C        -- Put in a trivial but scaled equation
         DO 180 J = 1, 7
            A(ISE,J) = ZERO
  180    CONTINUE
         A(ISE,O) = SCALE
         IF (BBNDX.EQ.ZERO .AND. BBNDY.EQ.ZERO) THEN
            USOL = HALF*(CBNDX/ABNDX+CBNDY/ABNDY)
         ELSE IF (BBNDX.EQ.ZERO) THEN
            USOL = CBNDX/ABNDX
         ELSE IF (BBNDY.EQ.ZERO) THEN
            USOL = CBNDY/ABNDY
         END IF
         RHS(ISE) = USOL*SCALE
         NEUMAN = .FALSE.
      ELSE
C        -- Both Boundary Conditions are Mixed
         A(ISE,W) = A(ISE,W) + A(ISE,E)
         A(ISE,N) = A(ISE,N) + A(ISE,S)
         A(ISE,O) = A(ISE,O) - TWO*HX*ABNDX*A(ISE,E)/BBNDX -
     *              TWO*HY*ABNDY*A(ISE,S)/BBNDY
         RHS(ISE) = RHS(ISE) - TWO*HX*CBNDX*A(ISE,E)/BBNDX -
     *              TWO*HY*CBNDY*A(ISE,S)/BBNDY
         A(ISE,S) = ZERO
         A(ISE,E) = ZERO
         IF (ABNDX.NE.ZERO) NEUMAN = .FALSE.
         IF (ABNDY.NE.ZERO) NEUMAN = .FALSE.
         IF (A(ISE,SE).NE.ZERO .OR. A(ISE,NW).NE.ZERO) THEN
            IF (IFAIL1.NE.1) THEN
               REC = ' ** Mixed Derivative in Equation '//
     *               'and Derivative in Boundary Condition'
               CALL X04BAF(NERR,REC)
               WRITE (REC,FMT='(A,1P,2E12.5)')
     *           ' ** South-East Corner, (x,y) = ', X, Y
               CALL X04BAF(NERR,REC)
            END IF
            IERROR = 2
            RETURN
         END IF
      END IF
C
C     Bottom Left Corner (SW) --
C
      X = XMIN
      Y = YMIN
      CALL BNDY(X,Y,ABNDX,BBNDX,CBNDX,BOTTOM)
      IF (ABNDX.EQ.ZERO .AND. BBNDX.EQ.ZERO) THEN
         IF (IFAIL1.NE.1) THEN
            CALL X04BAF(NERR,' ** Null Boundary Condition')
            WRITE (REC,FMT='(A,1P,2E12.5)')
     *        ' ** Bottom Boundary, Left End, (x,y) = ', X, Y
            CALL X04BAF(NERR,REC)
         END IF
         IERROR = 3
         RETURN
      END IF
      CALL BNDY(X,Y,ABNDY,BBNDY,CBNDY,LEFT)
      IF (ABNDY.EQ.ZERO .AND. BBNDY.EQ.ZERO) THEN
         IF (IFAIL1.NE.1) THEN
            CALL X04BAF(NERR,' ** Null Boundary Condition')
            WRITE (REC,FMT='(A,1P,2E12.5)')
     *        ' ** Left Boundary, Bottom End, (x,y) = ', X, Y
            CALL X04BAF(NERR,REC)
         END IF
         IERROR = 3
         RETURN
      END IF
      IF (BBNDX.EQ.ZERO .OR. BBNDY.EQ.ZERO) THEN
C        -- One of the Boundary Conditions is Dirichlet
C        -- Put in a trivial but scaled equation
         DO 200 J = 1, 7
            A(ISW,J) = ZERO
  200    CONTINUE
         A(ISW,O) = SCALE
         IF (BBNDX.EQ.ZERO .AND. BBNDY.EQ.ZERO) THEN
            USOL = HALF*(CBNDX/ABNDX+CBNDY/ABNDY)
         ELSE IF (BBNDX.EQ.ZERO) THEN
            USOL = CBNDX/ABNDX
         ELSE IF (BBNDY.EQ.ZERO) THEN
            USOL = CBNDY/ABNDY
         END IF
         RHS(ISW) = USOL*SCALE
         NEUMAN = .FALSE.
      ELSE
C        -- Both Boundary Conditions are Mixed
         A(ISW,E) = A(ISW,E) + A(ISW,W)
         A(ISW,N) = A(ISW,N) + A(ISW,S)
         A(ISW,O) = A(ISW,O) - TWO*HX*ABNDX*A(ISW,W)/BBNDX -
     *              TWO*HY*ABNDY*A(ISW,S)/BBNDY
         RHS(ISW) = RHS(ISW) - TWO*HX*CBNDX*A(ISW,W)/BBNDX -
     *              TWO*HY*CBNDY*A(ISW,S)/BBNDY
         A(ISW,S) = ZERO
         A(ISW,W) = ZERO
         IF (ABNDX.NE.ZERO) NEUMAN = .FALSE.
         IF (ABNDY.NE.ZERO) NEUMAN = .FALSE.
         IF (A(ISW,SE).NE.ZERO .OR. A(ISW,NW).NE.ZERO) THEN
            IF (IFAIL1.NE.1) THEN
               REC = ' ** Mixed Derivative in Equation '//
     *               'and Derivative in Boundary Condition'
               CALL X04BAF(NERR,REC)
               WRITE (REC,FMT='(A,1P,2E12.5)')
     *           ' ** South-West Corner, (x,y) = ', X, Y
               CALL X04BAF(NERR,REC)
            END IF
            IERROR = 2
            RETURN
         END IF
      END IF
C
C     Top Right Corner (NE) --
C
      X = XMAX
      Y = YMAX
      CALL BNDY(X,Y,ABNDX,BBNDX,CBNDX,TOP)
      IF (ABNDX.EQ.ZERO .AND. BBNDX.EQ.ZERO) THEN
         IF (IFAIL1.NE.1) THEN
            CALL X04BAF(NERR,' ** Null Boundary Condition')
            WRITE (REC,FMT='(A,1P,2E12.5)')
     *        ' ** Top Boundary, Right End, (x,y) = ', X, Y
            CALL X04BAF(NERR,REC)
         END IF
         IERROR = 3
         RETURN
      END IF
      CALL BNDY(X,Y,ABNDY,BBNDY,CBNDY,RIGHT)
      IF (ABNDY.EQ.ZERO .AND. BBNDY.EQ.ZERO) THEN
         IF (IFAIL1.NE.1) THEN
            CALL X04BAF(NERR,' ** Null Boundary Condition')
            WRITE (REC,FMT='(A,1P,2E12.5)')
     *        ' ** Right Boundary, Top End, (x,y) = ', X, Y
            CALL X04BAF(NERR,REC)
         END IF
         IERROR = 3
         RETURN
      END IF
      IF (BBNDX.EQ.ZERO .OR. BBNDY.EQ.ZERO) THEN
C        -- One of the Boundary Conditions is Dirichlet
C        -- Put in a trivial but scaled equation
         DO 220 J = 1, 7
            A(INE,J) = ZERO
  220    CONTINUE
         A(INE,O) = SCALE
         IF (BBNDX.EQ.ZERO .AND. BBNDY.EQ.ZERO) THEN
            USOL = HALF*(CBNDX/ABNDX+CBNDY/ABNDY)
         ELSE IF (BBNDX.EQ.ZERO) THEN
            USOL = CBNDX/ABNDX
         ELSE IF (BBNDY.EQ.ZERO) THEN
            USOL = CBNDY/ABNDY
         END IF
         RHS(INE) = USOL*SCALE
         NEUMAN = .FALSE.
      ELSE
C        -- Both Boundary Conditions are Mixed
         A(INE,W) = A(INE,W) + A(INE,E)
         A(INE,S) = A(INE,S) + A(INE,N)
         A(INE,O) = A(INE,O) - TWO*HX*ABNDX*A(INE,E)/BBNDX -
     *              TWO*HY*ABNDY*A(INE,N)/BBNDY
         RHS(INE) = RHS(INE) - TWO*HX*CBNDX*A(INE,E)/BBNDX -
     *              TWO*HY*CBNDY*A(INE,N)/BBNDY
         A(INE,N) = ZERO
         A(INE,E) = ZERO
         IF (ABNDX.NE.ZERO) NEUMAN = .FALSE.
         IF (ABNDY.NE.ZERO) NEUMAN = .FALSE.
         IF (A(INE,SE).NE.ZERO .OR. A(INE,NW).NE.ZERO) THEN
            IF (IFAIL1.NE.1) THEN
               REC = ' ** Mixed Derivative in Equation '//
     *               'and Derivative in Boundary Condition'
               CALL X04BAF(NERR,REC)
               WRITE (REC,FMT='(A,1P,2E12.5)')
     *           ' ** North-East Corner, (x,y) = ', X, Y
               CALL X04BAF(NERR,REC)
            END IF
            IERROR = 2
            RETURN
         END IF
      END IF
C
C     Top Left Corner(NW) --
C
      X = XMIN
      Y = YMAX
      CALL BNDY(X,Y,ABNDX,BBNDX,CBNDX,TOP)
      IF (ABNDX.EQ.ZERO .AND. BBNDX.EQ.ZERO) THEN
         IF (IFAIL1.NE.1) THEN
            CALL X04BAF(NERR,' ** Null Boundary Condition')
            WRITE (REC,FMT='(A,1P,2E12.5)')
     *        ' ** Top Boundary, Left End, (x,y) = ', X, Y
            CALL X04BAF(NERR,REC)
         END IF
         IERROR = 3
         RETURN
      END IF
      CALL BNDY(X,Y,ABNDY,BBNDY,CBNDY,LEFT)
      IF (ABNDY.EQ.ZERO .AND. BBNDY.EQ.ZERO) THEN
         IF (IFAIL1.NE.1) THEN
            CALL X04BAF(NERR,' ** Null Boundary Condition')
            WRITE (REC,FMT='(A,1P,2E12.5)')
     *        ' ** Left Boundary, Top End, (x,y) = ', X, Y
            CALL X04BAF(NERR,REC)
         END IF
         IERROR = 3
         RETURN
      END IF
      IF (BBNDX.EQ.ZERO .OR. BBNDY.EQ.ZERO) THEN
C        -- One of the Boundary Conditions is Dirichlet
C        -- Put in a trivial but scaled equation
         DO 240 J = 1, 7
            A(INW,J) = ZERO
  240    CONTINUE
         A(INW,O) = SCALE
         IF (BBNDX.EQ.ZERO .AND. BBNDY.EQ.ZERO) THEN
            USOL = HALF*(CBNDX/ABNDX+CBNDY/ABNDY)
         ELSE IF (BBNDX.EQ.ZERO) THEN
            USOL = CBNDX/ABNDX
         ELSE IF (BBNDY.EQ.ZERO) THEN
            USOL = CBNDY/ABNDY
         END IF
         RHS(INW) = USOL*SCALE
         NEUMAN = .FALSE.
      ELSE
C        -- Both Boundary Conditions are Mixed
         A(INW,E) = A(INW,E) + A(INW,W)
         A(INW,S) = A(INW,S) + A(INW,N)
         A(INW,O) = A(INW,O) - TWO*HX*ABNDX*A(INW,W)/BBNDX -
     *              TWO*HY*ABNDY*A(INW,N)/BBNDY
         RHS(INW) = RHS(INW) - TWO*HX*CBNDX*A(INW,W)/BBNDX -
     *              TWO*HY*CBNDY*A(INW,N)/BBNDY
         A(INW,N) = ZERO
         A(INW,W) = ZERO
         IF (ABNDX.NE.ZERO) NEUMAN = .FALSE.
         IF (ABNDY.NE.ZERO) NEUMAN = .FALSE.
         IF (A(INW,SE).NE.ZERO .OR. A(INW,NW).NE.ZERO) THEN
            IF (IFAIL1.NE.1) THEN
               REC = ' ** Mixed Derivative in Equation '//
     *               'and Derivative in Boundary Condition'
               CALL X04BAF(NERR,REC)
               WRITE (REC,FMT='(A,1P,2E12.5)')
     *           ' ** North-West Corner, (x,y) = ', X, Y
               CALL X04BAF(NERR,REC)
            END IF
            IERROR = 2
            RETURN
         END IF
      END IF
C
C     Check the 7-diagonal equations for diagonal dominance --
C
      FUDGE = 1.0D0 + 10.D0*EPS
      DIADOM = .TRUE.
      DO 260 I = 1, NGX*NGY
         SUM = ABS(A(I,S)) + ABS(A(I,SE)) + ABS(A(I,W)) + ABS(A(I,E)) +
     *         ABS(A(I,NW)) + ABS(A(I,N))
         IF (SUM.GT.FUDGE*ABS(A(I,O))) DIADOM = .FALSE.
  260 CONTINUE
      IF ( .NOT. DIADOM) IERROR = 6
      IF (NEUMAN) IERROR = 5
      IF ( .NOT. ELLPTC) IERROR = 4
C
      RETURN
      END
