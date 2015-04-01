      SUBROUTINE E01SAW(NST,PX,PY,X,Y,IADJ,IEND,I1,I2,I3)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 13B REVISED. IER-656 (AUG 1988).
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Original name: TRFIND.
C     This routine locates a point P in a Thiessen triangu-
C     lation, returning the vertex indices of a triangle which
C     contains P.  E01SAW is part of an interpolation package
C     which provides subroutines for creating the mesh.
C
C     Input Parameters -    NST - index of node at which E01SAW
C                             begins search.  Search time
C                             depends on the proximity of
C                             NST to P.
C
C                     PX,PY - X and Y-coordinates of the
C                             point to be located.
C
C                       X,Y - vectors of coordinates of
C                             nodes in the mesh.  (X(I),Y(I))
C                             defines node I for I = 1,...,N
C                             where N .ge. 3.
C
C                      IADJ - set of adjacency lists of
C                             nodes in the mesh.
C
C                      IEND - pointers to the ends of
C                             adjacency lists in IADJ for
C                             each node in the mesh.
C
C     IADJ and IEND may be created by E01SAY.
C
C     Input parameters are not altered by this routine.
C
C     Output Parameters - I1,I2,I3 - vertex indices in counter-
C                                clockwise order - vertices
C                                of a triangle containing P
C                                if P is an interior node.
C                                If P is outside of the
C                                boundary of the mesh, I1
C                                and I2 are the first (right
C                                -most) and last (leftmost)
C                                nodes which are visible
C                                from P, and I3 = 0.  If P
C                                is on the mesh boundary, I1
C                                and I2 are the first (right
C                                -most) and last (leftmost)
C                                nodes which are visible
C                                from P, and I3 = -1. If P
C                                and all of the nodes lie on
C                                a single line then I1 = I2
C                                = I3 = 0.
C
C     ***********************************************************
C
C     Local Parameters -
C
C     XP,YP =     local variables containing PX and PY
C     N0,N1,N2 =  nodes in counterclockwise order defining a
C               cone (with vertex N0) containing P
C     N3,N4 =     nodes opposite N1-N2 and N2-N1, respectively
C     INDX,IND =  indices for IADJ
C     NF,NL =     first and last neighbors of N0 in IADJ, or
C               first (rightmost) and last (leftmost) nodes
C               visible from P when P is outside the
C               boundary
C     NEXT =      candidate for I1 or I2 when P is outside of
C               the boundary
C     LEFT =      statement function which computes the sign of
C               a cross product (Z-component).  LEFT(X1,...,
C               Y0) = .TRUE. iff (X0,Y0) is on or to the
C               left of the vector from (X1,Y1) to (X2,Y2).
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  PX, PY
      INTEGER           I1, I2, I3, NST
C     .. Array Arguments ..
      DOUBLE PRECISION  X(*), Y(*)
      INTEGER           IADJ(*), IEND(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  X0, X1, X2, XP, Y0, Y1, Y2, YP
      INTEGER           IND, INDX, N0, N1, N2, N3, N4, NEXT, NF, NL
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Statement Functions ..
      LOGICAL           LEFT
C     .. Statement Function definitions ..
      LEFT(X1,Y1,X2,Y2,X0,Y0) = (X2-X1)*(Y0-Y1) .GE. (X0-X1)*(Y2-Y1)
C     .. Executable Statements ..
      XP = PX
      YP = PY
C
C     Initialize variables and find a cone containing P
C
      N0 = MAX(NST,1)
   20 CONTINUE
      INDX = IEND(N0)
      NL = IADJ(INDX)
      INDX = 1
      IF (N0.NE.1) INDX = IEND(N0-1) + 1
      NF = IADJ(INDX)
      N1 = NF
      IF (NL.NE.0) THEN
   40    CONTINUE
C
C        N0 is an interior node.  Find N1.
C
         IF ( .NOT. LEFT(X(N0),Y(N0),X(N1),Y(N1),XP,YP)) THEN
            INDX = INDX + 1
            N1 = IADJ(INDX)
            IF (N1.EQ.NL) THEN
               GO TO 160
            ELSE
               GO TO 40
            END IF
         END IF
      ELSE
C
C        N0 is a boundary node.  Set NL to the last nonzero
C        neighbor of N0.
C
         IND = IEND(N0) - 1
         NL = IADJ(IND)
         IF ( .NOT. LEFT(X(N0),Y(N0),X(NF),Y(NF),XP,YP)) THEN
            GO TO 120
         ELSE IF ( .NOT. LEFT(X(NL),Y(NL),X(N0),Y(N0),XP,YP)) THEN
            GO TO 100
         END IF
      END IF
   60 CONTINUE
C
C     P is to the left of arc N0-N1.  Initialize N2 to the next
C     neighbor of N0.
C
      INDX = INDX + 1
      N2 = IADJ(INDX)
      IF (LEFT(X(N0),Y(N0),X(N2),Y(N2),XP,YP)) THEN
         N1 = N2
         IF (N1.NE.NL) GO TO 60
      ELSE
         GO TO 180
      END IF
      IF (LEFT(X(N0),Y(N0),X(NF),Y(NF),XP,YP)) THEN
         IF (XP.NE.X(N0) .OR. YP.NE.Y(N0)) THEN
   80       CONTINUE
C
C           P is left of or on arcs N0-NB for all neighbors NB
C           of N0.
C           All points are collinear iff P is left of NB-N0 for
C           all neighbors NB of N0.  Search the neighbors of N0
C           in reverse order.  Note -- N1 = NL and INDX points to
C           NL.
C
            IF (LEFT(X(N1),Y(N1),X(N0),Y(N0),XP,YP)) THEN
               IF (N1.EQ.NF) THEN
                  GO TO 140
               ELSE
                  INDX = INDX - 1
                  N1 = IADJ(INDX)
                  GO TO 80
               END IF
            END IF
         END IF
C
C        P is to the right of N1-N0, or P=N0.  Set N0 to N1 and
C        start over.
C
         N0 = N1
         GO TO 20
      ELSE
         GO TO 160
      END IF
C
C     P is outside the boundary and N0 is the rightmost
C     visible boundary node
C
  100 I1 = N0
      GO TO 300
C
C     P is outside the boundary
C
  120 NL = N0
      GO TO 280
C
C     All points are collinear
C
  140 I1 = 0
      I2 = 0
      I3 = 0
      RETURN
C
C     P is between arcs N0-N1 and N0-NF
C
  160 N2 = NF
C
C     P is contained in a cone defined by line segments N0-N1
C     and N0-N2 where N1 is adjacent to N2
C
  180 N3 = N0
  200 CONTINUE
      IF (LEFT(X(N1),Y(N1),X(N2),Y(N2),XP,YP)) THEN
C        The algorithm believes that P is contained in
C        triangle (N1,N2,N3). However, if N1, N2, N3 and P
C        are almost collinear then P may in fact be far
C        outside the triangle, the algorithm being fooled
C        by rounding errors. The following code tests
C        whether P lies inside a rectangle containing the
C        triangle. If not, the triangle is rejected, and
C        the algorithm continues.
         IF (XP.GE.MIN(X(N1),X(N2),X(N3)) .AND.
     *       XP.LE.MAX(X(N1),X(N2),X(N3)) .AND.
     *       YP.GE.MIN(Y(N1),Y(N2),Y(N3)) .AND.
     *       YP.LE.MAX(Y(N1),Y(N2),Y(N3))) GO TO 260
      END IF
C
C     Set N4 to the first neighbor of N2 following N1
C
      INDX = IEND(N2)
      IF (IADJ(INDX).NE.N1) THEN
  220    CONTINUE
C
C        N1 is not the last neighbor of N2
C
         INDX = INDX - 1
         IF (IADJ(INDX).NE.N1) GO TO 220
         N4 = IADJ(INDX+1)
         IF (N4.EQ.0) GO TO 240
      ELSE
C
C        N1 is the last neighbor of N2.
C        Set N4 to the first neighbor.
C
         INDX = 1
         IF (N2.NE.1) INDX = IEND(N2-1) + 1
         N4 = IADJ(INDX)
      END IF
C
C     Define a new arc N1-N2 which intersects the line
C     segment N0-P
C
      IF (LEFT(X(N0),Y(N0),X(N4),Y(N4),XP,YP)) THEN
         N3 = N1
         N1 = N4
      ELSE
         N3 = N2
         N2 = N4
      END IF
      GO TO 200
C
C     P is outside the boundary
C
  240 NF = N2
      NL = N1
      GO TO 280
C
C     P is in the triangle (N1,N2,N3) and not on N2-N3.  If
C     N3-N1 or N1-N2 is a boundary arc containing P, treat P
C     as exterior.
C
  260 INDX = IEND(N1)
      IF (IADJ(INDX).EQ.0) THEN
C
C        N1 is a boundary node.  N3-N1 is a boundary arc iff N3
C        is the last nonzero neighbor of N1.
C
         IF (N3.EQ.IADJ(INDX-1)) THEN
C
C           N3-N1 is a boundary arc
C
            IF (LEFT(X(N1),Y(N1),X(N3),Y(N3),XP,YP)) THEN
C
C              P lies on N1-N3
C
               I1 = N1
               I2 = N3
               I3 = -1
               RETURN
            END IF
         END IF
C
C        N3-N1 is not a boundary arc containing P.  N1-N2 is a
C        boundary arc iff N2 is the first neighbor of N1.
C
         INDX = 1
         IF (N1.NE.1) INDX = IEND(N1-1) + 1
         IF (N2.EQ.IADJ(INDX)) THEN
C
C           N1-N2 is a boundary arc
C
            IF (LEFT(X(N2),Y(N2),X(N1),Y(N1),XP,YP)) THEN
C
C              P lies on N1-N2
C
               I1 = N2
               I2 = N1
               I3 = -1
               RETURN
            END IF
         END IF
      END IF
C
C     P does not lie on a boundary arc.
C
      I1 = N1
      I2 = N2
      I3 = N3
      RETURN
  280 CONTINUE
C
C     NF and NL are adjacent boundary nodes which are visible
C     from P.  Find the first visible boundary node.
C     Set next to the first neighbor of NF.
C
      INDX = 1
      IF (NF.NE.1) INDX = IEND(NF-1) + 1
      NEXT = IADJ(INDX)
      IF ( .NOT. LEFT(X(NF),Y(NF),X(NEXT),Y(NEXT),XP,YP)) THEN
         NF = NEXT
         GO TO 280
      END IF
C
C     NF is the first (rightmost) visible boundary node
C
      I1 = NF
  300 CONTINUE
C
C     Find the last visible boundary node.  NL is the first
C     candidate for I2.
C     Set next to the last neighbor of NL.
C
      INDX = IEND(NL) - 1
      NEXT = IADJ(INDX)
      IF ( .NOT. LEFT(X(NEXT),Y(NEXT),X(NL),Y(NL),XP,YP)) THEN
         NL = NEXT
         GO TO 300
      END IF
C
C     NL is the last (leftmost) visible boundary node
C
      I2 = NL
      I3 = 0
      RETURN
      END
