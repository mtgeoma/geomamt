      SUBROUTINE E01SBY(X,Y,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,ZX1,ZX2,ZX3,ZY1,
     *                  ZY2,ZY3,IFLAG,W,WX,WY,IER)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Given function values and first partial derivatives at
C     the three vertices of a triangle, this routine determines
C     a function W which agrees with the given data, returning
C     the value and (optionally) first partial derivatives of W
C     at a point (X,Y) in the triangle.  The interpolation
C     method is exact for quadratic polynomial data.  The
C     triangle is partitioned into three subtriangles with
C     equal areas.  W is cubic in each subtriangle and along
C     the edges, but has only one continuous derivative across
C     edges.  The normal derivative of W varies linearly along
C     each outer edge.  The values and partial derivatives of W
C     along a triangle edge depend only on the data values at
C     the endpoints of the edge.  Thus the method yields C-1
C     continuity when used to interpolate over a triangular
C     grid.  This algorithm is due to C. L. Lawson.
C
C     Input Parameters -   X,Y - coordinates of a point at which
C                            W is to be evaluated.
C
C        X1,X2,X3,Y1,Y2,Y3 - coordinates of the vertices of
C                            a triangle containing (X,Y).
C
C                 Z1,Z2,Z3 - function values at the vertices
C                            to be interpolated.
C
C              ZX1,ZX2,ZX3 - X-derivative values at the
C                            vertices.
C
C              ZY1,ZY2,ZY3 - Y-derivative values at the
C                            vertices.
C
C                    IFLAG - option indicator
C                            IFLAG = 0 if only W is to be
C                                      computed.
C                            IFLAG = 1 if W, WX, and WY are
C                                      to be returned.
C
C     Input parameters are not altered by this routine.
C
C     Output Parameters -   W - estimated value of the interp-
C                           olatory function at (X,Y) if
C                           IER = 0.  Otherwise W = 0.
C
C                   WX,WY - partial derivatives of W at
C                           (X,Y) if IER = 0 and IFLAG = 1,
C                           unchanged if IFLAG .ne. 1, zero
C                           if IER .ne. 0 and IFLAG = 1.
C
C                     IER - error indicator
C                           IER = 0 if no errors were
C                                   encountered.
C                           IER = 1 if the vertices of the
C                                   triangle are collinear.
C
C     ***********************************************************
C
C     Local Parameters -
C
C     I =               do-loop index
C     IP1,IP2,IP3 =     permuted indices for computing RO, ROX,
C                       and ROY
C     U(K) =            X-component of the vector representing
C                       the side opposite vertex K
C     V(K) =            Y-component of the vector representing
C                       the side opposite vertex K
C     SL(K) =           square of the length of the side
C                       opposite vertex K
C     AREA =            twice the area of the triangle
C     XP,YP =           X-X1, Y-Y1
C     R(K) =            K-th barycentric coordinate
C     RX(K),RY(K) =     X,Y partial derivatives of R(K)
C     PHI(K)            R(K-1)*R(K+1) -- quadratic
C     PHIX(K),PHIY(K) = X,Y partials of PHI(K)
C     RMIN =            min(R1,R2,R3)
C     C1,C2 =           factors for computing RO
C     RO(K) =           factors for computing G -- cubic
C                       correction terms
C     ROX(K),ROY(K) =   X,Y partials of RO(K)
C     F(K) =            factors for computing G, GX, and GY --
C                       constant
C     G(K) =            factors for computing the cardinal
C                       functions -- cubic
C     GX(K),GY(K) =     X,Y partials of G(K)
C     P(K) =            G(K) + PHI(K)
C     PX(K),PY(K) =     X,Y partials of P(K)
C     Q(K) =            G(K) - PHI(K)
C     QX(K),QY(K) =     X,Y partials of Q(K)
C     A(K) =            cardinal function whose coefficient is
C                       Z(K)
C     AX(K),AY(K) =     X,Y partials of A(K) -- cardinal
C                       functions for WX and WY
C     B(K) =            twice the cardinal function whose
C                       coefficient is ZX(K)
C     BX(K),BY(K) =     X,Y partials of B(K)
C     C(K) =            twice the cardinal function whose
C                       coefficient is ZY(K)
C     CX(K),CY(K) =     X,Y partials of C(K)
C
C     .. Parameters ..
      DOUBLE PRECISION  FIVE, THREE, TWO, ZERO
      PARAMETER         (FIVE=5.0D0,THREE=3.0D0,TWO=2.0D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  W, WX, WY, X, X1, X2, X3, Y, Y1, Y2, Y3, Z1, Z2,
     *                  Z3, ZX1, ZX2, ZX3, ZY1, ZY2, ZY3
      INTEGER           IER, IFLAG
C     .. Local Scalars ..
      DOUBLE PRECISION  AREA, C1, C2, RMIN, XP, YP
      INTEGER           I, IP1, IP2, IP3
C     .. Local Arrays ..
      DOUBLE PRECISION  A(3), AX(3), AY(3), B(3), BX(3), BY(3), C(3),
     *                  CX(3), CY(3), F(3), G(3), GX(3), GY(3), P(3),
     *                  PHI(3), PHIX(3), PHIY(3), PX(3), PY(3), Q(3),
     *                  QX(3), QY(3), R(3), RO(3), ROX(3), ROY(3),
     *                  RX(3), RY(3), SL(3), U(3), V(3)
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
      U(1) = X3 - X2
      U(2) = X1 - X3
      U(3) = X2 - X1
C
      V(1) = Y3 - Y2
      V(2) = Y1 - Y3
      V(3) = Y2 - Y1
C
      DO 20 I = 1, 3
         SL(I) = U(I)*U(I) + V(I)*V(I)
   20 CONTINUE
C
C     Area = 3-1 x 3-2
C
      AREA = U(1)*V(2) - U(2)*V(1)
      IF (AREA.EQ.ZERO) THEN
C
C        Vertices are collinear
C
         IER = 1
         W = ZERO
         IF (IFLAG.EQ.1) THEN
            WX = ZERO
            WY = ZERO
         END IF
      ELSE
C
C        R(1) = (2-3 X 2-(X,Y))/AREA, R(2) = (1-(X,Y) X 1-3)/AREA,
C        R(3) = (1-2 X 1-(X,Y))/AREA
C
         R(1) = (U(1)*(Y-Y2)-V(1)*(X-X2))/AREA
         XP = X - X1
         YP = Y - Y1
         R(2) = (U(2)*YP-V(2)*XP)/AREA
         R(3) = (U(3)*YP-V(3)*XP)/AREA
         IER = 0
C
         PHI(1) = R(2)*R(3)
         PHI(2) = R(3)*R(1)
         PHI(3) = R(1)*R(2)
C
         RMIN = MIN(R(1),R(2),R(3))
         IF (RMIN.EQ.R(1)) THEN
            IP1 = 1
            IP2 = 2
            IP3 = 3
         ELSE IF (RMIN.NE.R(2)) THEN
            IP1 = 3
            IP2 = 1
            IP3 = 2
         ELSE
            IP1 = 2
            IP2 = 3
            IP3 = 1
         END IF
C
         C1 = RMIN*RMIN/TWO
         C2 = RMIN/THREE
         RO(IP1) = (PHI(IP1)+FIVE*C1/THREE)*R(IP1) - C1
         RO(IP2) = C1*(R(IP3)-C2)
         RO(IP3) = C1*(R(IP2)-C2)
C
         F(1) = THREE*(SL(2)-SL(3))/SL(1)
         F(2) = THREE*(SL(3)-SL(1))/SL(2)
         F(3) = THREE*(SL(1)-SL(2))/SL(3)
C
         G(1) = (R(2)-R(3))*PHI(1) + F(1)*RO(1) - RO(2) + RO(3)
         G(2) = (R(3)-R(1))*PHI(2) + F(2)*RO(2) - RO(3) + RO(1)
         G(3) = (R(1)-R(2))*PHI(3) + F(3)*RO(3) - RO(1) + RO(2)
C
         DO 40 I = 1, 3
            P(I) = G(I) + PHI(I)
            Q(I) = G(I) - PHI(I)
   40    CONTINUE
C
         A(1) = R(1) + G(3) - G(2)
         A(2) = R(2) + G(1) - G(3)
         A(3) = R(3) + G(2) - G(1)
C
         B(1) = U(3)*P(3) + U(2)*Q(2)
         B(2) = U(1)*P(1) + U(3)*Q(3)
         B(3) = U(2)*P(2) + U(1)*Q(1)
C
         C(1) = V(3)*P(3) + V(2)*Q(2)
         C(2) = V(1)*P(1) + V(3)*Q(3)
         C(3) = V(2)*P(2) + V(1)*Q(1)
C
C        W is a linear combination of the cardinal functions
C
         W = A(1)*Z1 + A(2)*Z2 + A(3)*Z3 + (B(1)*ZX1+B(2)*ZX2+B(3)
     *       *ZX3+C(1)*ZY1+C(2)*ZY2+C(3)*ZY3)/TWO
         IF (IFLAG.EQ.1) THEN
C
C           Compute WX and WY
C
            DO 60 I = 1, 3
               RX(I) = -V(I)/AREA
               RY(I) = U(I)/AREA
   60       CONTINUE
            PHIX(1) = R(2)*RX(3) + RX(2)*R(3)
            PHIY(1) = R(2)*RY(3) + RY(2)*R(3)
            PHIX(2) = R(3)*RX(1) + RX(3)*R(1)
            PHIY(2) = R(3)*RY(1) + RY(3)*R(1)
            PHIX(3) = R(1)*RX(2) + RX(1)*R(2)
            PHIY(3) = R(1)*RY(2) + RY(1)*R(2)
C
            ROX(IP1) = RX(IP1)*(PHI(IP1)+FIVE*C1) + R(IP1)*(PHIX(IP1)
     *                 -RX(IP1))
            ROY(IP1) = RY(IP1)*(PHI(IP1)+FIVE*C1) + R(IP1)*(PHIY(IP1)
     *                 -RY(IP1))
            ROX(IP2) = RX(IP1)*(PHI(IP2)-C1) + C1*RX(IP3)
            ROY(IP2) = RY(IP1)*(PHI(IP2)-C1) + C1*RY(IP3)
            ROX(IP3) = RX(IP1)*(PHI(IP3)-C1) + C1*RX(IP2)
            ROY(IP3) = RY(IP1)*(PHI(IP3)-C1) + C1*RY(IP2)
C
            GX(1) = (RX(2)-RX(3))*PHI(1) + (R(2)-R(3))*PHIX(1) + F(1)
     *              *ROX(1) - ROX(2) + ROX(3)
            GY(1) = (RY(2)-RY(3))*PHI(1) + (R(2)-R(3))*PHIY(1) + F(1)
     *              *ROY(1) - ROY(2) + ROY(3)
            GX(2) = (RX(3)-RX(1))*PHI(2) + (R(3)-R(1))*PHIX(2) + F(2)
     *              *ROX(2) - ROX(3) + ROX(1)
            GY(2) = (RY(3)-RY(1))*PHI(2) + (R(3)-R(1))*PHIY(2) + F(2)
     *              *ROY(2) - ROY(3) + ROY(1)
            GX(3) = (RX(1)-RX(2))*PHI(3) + (R(1)-R(2))*PHIX(3) + F(3)
     *              *ROX(3) - ROX(1) + ROX(2)
            GY(3) = (RY(1)-RY(2))*PHI(3) + (R(1)-R(2))*PHIY(3) + F(3)
     *              *ROY(3) - ROY(1) + ROY(2)
C
            DO 80 I = 1, 3
               PX(I) = GX(I) + PHIX(I)
               PY(I) = GY(I) + PHIY(I)
               QX(I) = GX(I) - PHIX(I)
               QY(I) = GY(I) - PHIY(I)
   80       CONTINUE
C
            AX(1) = RX(1) + GX(3) - GX(2)
            AY(1) = RY(1) + GY(3) - GY(2)
            AX(2) = RX(2) + GX(1) - GX(3)
            AY(2) = RY(2) + GY(1) - GY(3)
            AX(3) = RX(3) + GX(2) - GX(1)
            AY(3) = RY(3) + GY(2) - GY(1)
C
            BX(1) = U(3)*PX(3) + U(2)*QX(2)
            BY(1) = U(3)*PY(3) + U(2)*QY(2)
            BX(2) = U(1)*PX(1) + U(3)*QX(3)
            BY(2) = U(1)*PY(1) + U(3)*QY(3)
            BX(3) = U(2)*PX(2) + U(1)*QX(1)
            BY(3) = U(2)*PY(2) + U(1)*QY(1)
C
            CX(1) = V(3)*PX(3) + V(2)*QX(2)
            CY(1) = V(3)*PY(3) + V(2)*QY(2)
            CX(2) = V(1)*PX(1) + V(3)*QX(3)
            CY(2) = V(1)*PY(1) + V(3)*QY(3)
            CX(3) = V(2)*PX(2) + V(1)*QX(1)
            CY(3) = V(2)*PY(2) + V(1)*QY(1)
C
C           WX and WY are linear combinations of the cardinal
C           functions
C
            WX = AX(1)*Z1 + AX(2)*Z2 + AX(3)*Z3 + (BX(1)*ZX1+BX(2)
     *           *ZX2+BX(3)*ZX3+CX(1)*ZY1+CX(2)*ZY2+CX(3)*ZY3)/TWO
            WY = AY(1)*Z1 + AY(2)*Z2 + AY(3)*Z3 + (BY(1)*ZX1+BY(2)
     *           *ZX2+BY(3)*ZX3+CY(1)*ZY1+CY(2)*ZY2+CY(3)*ZY3)/TWO
         END IF
      END IF
      RETURN
      END
