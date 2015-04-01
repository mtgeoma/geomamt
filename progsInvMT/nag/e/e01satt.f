      LOGICAL FUNCTION E01SAT(IN1,IN2,IO1,IO2,X,Y)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Original name: SWPTST.
C     This function decides whether or not to replace a
C     diagonal arc in a quadrilateral with the other diagonal.
C     The determination is based on the sizes of the angles
C     contained in the 2 triangles defined by the diagonal.
C     The diagonal is chosen to maximize the smallest of the
C     six angles over the two pairs of triangles.
C
C     Input Parameters -  IN1,IN2,IO1,IO2 - node indices of the
C                              four points defining the
C                              quadrilateral.  IO1 and IO2
C                              are currently connected by a
C                              diagonal arc.  This arc
C                              should be replaced by an arc
C                              connecting IN1, IN2 if the
C                              decision is made to swap.
C                              IN1,IO1,IO2 must be in
C                              counterclockwise order.
C
C                        X,Y - vectors of nodal coordinates.
C                              (X(I),Y(I)) are the coord-
C                              inates of node I for I = IN1,
C                              IN2, IO1, or IO2.
C
C     None of the input parameters are altered by this routine.
C
C     Output Parameter -  E01SAT - .TRUE. iff the arc connecting
C                              IO1 and IO2 is to be replaced
C
C     ***********************************************************
C
C     Local Parameters -
C
C     DX11,DY11 = X,Y coordinates of the vector IN1-IO1
C     DX12,DY12 = X,Y coordinates of the vector IN1-IO2
C     DX22,DY22 = X,Y coordinates of the vector IN2-IO2
C     DX21,DY21 = X,Y coordinates of the vector IN2-IO1
C     SIN1 =      cross product of the vectors IN1-IO1 and
C               IN1-IO2 -- proportional to sin(T1) where T1
C               is the angle at IN1 formed by the vectors
C     COS1 =      inner product of the vectors IN1-IO1 and
C               IN1-IO2 -- proportional to cos(T1)
C     SIN2 =      cross product of the vectors IN2-IO2 and
C               IN2-IO1 -- proportional to sin(T2) where T2
C               is the angle at IN2 formed by the vectors
C     COS2 =      inner product of the vectors IN2-IO2 and
C               IN2-IO1 -- proportional to cos(T2)
C     SIN12 =     SIN1*COS2 + COS1*SIN2 -- proportional to
C               sin(T1+T2)
C
C     .. Parameters ..
      DOUBLE PRECISION        ZERO
      PARAMETER               (ZERO=0.0D0)
C     .. Scalar Arguments ..
      INTEGER                 IN1, IN2, IO1, IO2
C     .. Array Arguments ..
      DOUBLE PRECISION        X(*), Y(*)
C     .. Local Scalars ..
      DOUBLE PRECISION        COS1, COS2, DX11, DX12, DX21, DX22, DY11,
     *                        DY12, DY21, DY22, SIN1, SIN12, SIN2
C     .. Executable Statements ..
      E01SAT = .FALSE.
C
C     Compute the vectors containing the angles T1, T2
C
      DX11 = X(IO1) - X(IN1)
      DX12 = X(IO2) - X(IN1)
      DX22 = X(IO2) - X(IN2)
      DX21 = X(IO1) - X(IN2)
C
      DY11 = Y(IO1) - Y(IN1)
      DY12 = Y(IO2) - Y(IN1)
      DY22 = Y(IO2) - Y(IN2)
      DY21 = Y(IO1) - Y(IN2)
C
C     Compute inner products
C
      COS1 = DX11*DX12 + DY11*DY12
      COS2 = DX22*DX21 + DY22*DY21
C
C     The diagonals should be swapped iff (T1+T2) .gt. 180
C     degrees.  The following two tests insure numerical
C     stability.
C
      IF (COS1.LT.ZERO .OR. COS2.LT.ZERO) THEN
         IF (COS1.GE.ZERO .OR. COS2.GE.ZERO) THEN
C
C           Compute vector cross products
C
            SIN1 = DX11*DY12 - DX12*DY11
            SIN2 = DX22*DY21 - DX21*DY22
            SIN12 = SIN1*COS2 + COS1*SIN2
            IF (SIN12.GE.ZERO) RETURN
         END IF
         E01SAT = .TRUE.
      END IF
      RETURN
      END
