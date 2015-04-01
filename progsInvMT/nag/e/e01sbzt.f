      SUBROUTINE E01SBZ(N,PX,PY,X,Y,Z,IADJ,IEND,ZXZY,IST,IFLAG,PZ,DZX,
     *                  DZY,IER)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14C REVISED. IER-875 (NOV 1990).
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Original name: INTRC1.
C     Adapted for NAG by H. Scullion Leics. Univ.
C
C     Given a triangulation of a set of points in the plane,
C     this routine determines a piecewise cubic function F(X,Y)
C     which interpolates a set of data values and partial
C     derivatives at the vertices.  F has continuous first
C     derivatives over the mesh and extends beyond the mesh
C     boundary allowing extrapolation.  Interpolation is exact
C     for quadratic data.  The value of F at (PX,PY) is
C     returned.  E01SBZ is part of an interpolation package
C     which provides routines to generate, update and plot the
C     mesh.
C
C     Input Parameters -     N - number of nodes in the mesh.
C                            N .ge. 3.
C
C                    PX,PY - coordinates of a point at which
C                            F is to be evaluated.
C
C                      X,Y - vectors of coordinates of the
C                            nodes in the mesh.
C
C                        Z - vector of data values at the
C                            nodes.
C
C                     IADJ - set of adjacency lists of nodes
C                            in the mesh.
C
C                     IEND - pointers to the ends of
C                            adjacency lists in IADJ for
C                            each node in the mesh.
C
C                     ZXZY - 2 by N array whose columns
C                            contain estimated partial der-
C                            ivatives at the nodes (X par-
C                            tials in the first row)
C
C                      IST - index of the starting node in
C                            the search for a triangle con-
C                            taining (PX,PY).  1 .le. IST
C                            .le. N.  The output value of
C                            IST from a previous call may
C                            be a good choice.
C
C     IADJ and IEND may be created by E01SAY and derivative
C     estimates are computed by E01SAZ.
C
C     Input parameters other than IST are not altered by this
C     routine.
C
C     Output Parameters - IST - index of one of the vertices of
C                           the triangle containing (PX,PY)
C                           unless IER .gt. 0.
C
C                      PZ - value of F at (PX,PY), or 0 if
C                           IER .gt. 0.
C
C                     IER - error indicator
C                           IER = 0 if no errors were
C                                   encountered.
C                           IER = 1 if N or IST is
C                                    out of range.
C                           IER = 2 if the nodes are col-
C                                    linear.
C                           IER = 3 if extrapolation was used.
C
C     ***********************************************************
C
C     Local Parameters -
C
C     NN =                      local copy of N
C     I1,I2,I3 =                vertices determined by E01SAW
C     IERR =                    error flag for calls to E01SBY
C     N1,N2 =                   endpoints of the closest bound-
C                             ary edge to P when P is out-
C                             side of the mesh boundary
C     INDX =                    IADJ index of N1 as a neighbor
C                             of N2
C     XP,YP =                   local copies of the coordinates
C                             of P=(PX,PY)
C     ZX1,ZY1,ZX2,ZY2,ZX3,ZY3 = X and Y derivatives at the
C                             vertices of a triangle T which
C                             contains P or at N1 and N2
C     X1,Y1,X2,Y2,X3,Y3 =       X,Y coordinates of the vertices
C                             of T or of N1 and N2
C     Z1,Z2,Z3 =                data values at the vertices of T
C     DP =                      inner product of N1-N2 and P-N2
C     U,V =                     X,Y coordinates of the vector
C                             N2-N1
C     XQ,YQ =                   X,Y coordinates of the closest
C                             boundary point to P when P is
C                             outside of the mesh boundary
C     R1,R2 =                   barycentric coordinates of Q
C                             with respect to the line seg-
C                             ment N2-N1 containing Q
C     A1,A2,B1,B2,C1,C2 =       cardinal functions for evaluat-
C                             ing the interpolatory surface
C                             at Q
C     F1,F2 =                   cubic factors used to compute
C                             the cardinal functions
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DZX, DZY, PX, PY, PZ
      INTEGER           IER, IFLAG, IST, N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N), Y(N), Z(N), ZXZY(2,N)
      INTEGER           IADJ(6*N), IEND(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, A2, B1, B2, C1, C2, DP, F1, F2, R1, R2, U,
     *                  V, X1, X2, X3, XP, XQ, Y1, Y2, Y3, YP, YQ, Z1,
     *                  Z2, Z3, ZX1, ZX2, ZX3, ZY1, ZY2, ZY3
      INTEGER           I1, I2, I3, IERR, INDX, N1, N2, NN
C     .. External Subroutines ..
      EXTERNAL          E01SAW, E01SBY
C     .. Executable Statements ..
      IER = 0
      NN = N
      PZ = ZERO
      DZX = ZERO
      DZY = ZERO
      IF (NN.LT.3) THEN
C
C        N or IST out of range.
C
         IER = 1
      ELSE
         IF (IST.LT.1 .OR. IST.GT.NN) IST = 1
         XP = PX
         YP = PY
C
C        Find a triangle containing P if P is within the mesh
C        boundary
C
         CALL E01SAW(IST,XP,YP,X,Y,IADJ,IEND,I1,I2,I3)
         IF (I1.NE.0) THEN
            IST = I1
            IF (I3.LE.0) THEN
               IF (I3.EQ.0) IER = 3
C
C              P is outside of the mesh boundary.  Extrapolate to P by
C              passing a linear function of one variable through the
C              value and directional derivative (in the direction
C              P-Q) of the interpolatory surface (E01SBY) at Q where
C              Q is the closest boundary point to P.
C
C              Determine Q by traversing the boundary starting from
C              the rightmost visible node I1.
C
               N2 = I1
   20          CONTINUE
C
C              Set N1 to the last nonzero neighbor of N2 and compute DP
C
               INDX = IEND(N2) - 1
               N1 = IADJ(INDX)
               X1 = X(N1)
               Y1 = Y(N1)
               X2 = X(N2)
               Y2 = Y(N2)
               DP = (X1-X2)*(XP-X2) + (Y1-Y2)*(YP-Y2)
               IF (DP.LE.ZERO) THEN
                  GO TO 40
               ELSE IF ((XP-X1)*(X2-X1)+(YP-Y1)*(Y2-Y1).LE.ZERO) THEN
                  N2 = N1
                  GO TO 20
               END IF
C
C              The closest boundary point Q lies on N2-N1.  Compute
C              partials at N1 and N2.
C
               ZX1 = ZXZY(1,N1)
               ZY1 = ZXZY(2,N1)
               ZX2 = ZXZY(1,N2)
               ZY2 = ZXZY(2,N2)
C
C              Compute Q, its barycentric coordinates, and the cardinal
C              functions for extrapolation
C
               U = X2 - X1
               V = Y2 - Y1
               R1 = DP/(U**2+V**2)
               R2 = ONE - R1
               XQ = R1*X1 + R2*X2
               YQ = R1*Y1 + R2*Y2
               F1 = R1*R1*R2
               F2 = R1*R2*R2
               A1 = R1 + (F1-F2)
               A2 = R2 - (F1-F2)
               B1 = U*F1
               B2 = -U*F2
               C1 = V*F1
               C2 = -V*F2
C
C              Compute the value of the interpolatory surface (E01SBY)
C              at Q
C
               PZ = A1*Z(N1) + A2*Z(N2) + B1*ZX1 + B2*ZX2 + C1*ZY1 +
     *              C2*ZY2
C
C              Compute the extrapolated value at P
C
               PZ = PZ + (R1*ZX1+R2*ZX2)*(XP-XQ) + (R1*ZY1+R2*ZY2)
     *              *(YP-YQ)
               RETURN
C
C              N2 is the closest boundary point to P.  Compute partial
C              derivatives at N2.
C
   40          CONTINUE
               ZX2 = ZXZY(1,N2)
               ZY2 = ZXZY(2,N2)
C
C              Compute extrapolated value at P
C
               PZ = Z(N2) + ZX2*(XP-X2) + ZY2*(YP-Y2)
               RETURN
            ELSE
C
C              Derivatives are user provided
C
               ZX1 = ZXZY(1,I1)
               ZX2 = ZXZY(1,I2)
               ZX3 = ZXZY(1,I3)
               ZY1 = ZXZY(2,I1)
               ZY2 = ZXZY(2,I2)
               ZY3 = ZXZY(2,I3)
C
C              Set local parameters for call to E01SBY
C
               X1 = X(I1)
               Y1 = Y(I1)
               X2 = X(I2)
               Y2 = Y(I2)
               X3 = X(I3)
               Y3 = Y(I3)
               Z1 = Z(I1)
               Z2 = Z(I2)
               Z3 = Z(I3)
               CALL E01SBY(XP,YP,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,ZX1,ZX2,ZX3,
     *                     ZY1,ZY2,ZY3,IFLAG,PZ,DZX,DZY,IERR)
               IF (IERR.EQ.0) RETURN
            END IF
         END IF
C
C        Nodes are collinear
C
         IER = 2
      END IF
      RETURN
      END
