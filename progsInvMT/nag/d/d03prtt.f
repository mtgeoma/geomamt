      SUBROUTINE D03PRT(XP,UP,IPTS,X,U,NPTS,NPDE,IFAIL,C)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C ----------------------------------------------------------------------
C     SPRINT cubic spline interpolation routine .
C     INTCUB routine from SPRINT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C       Inputs
C       ******
C       NPDE           : The number of PDEs (first dimension of the
C                        arrays U and UP ), assumed greater than 0.
C       NPTS           : Number of data points assumed greater than two
C       (X(I), U(K,I)) : Abscissae and ordinates of the data points.
C                        X is assumed to be strictly increasing.
C              K = 1,NPDE AND I = 1,NPTS.
C       XP(J)  J =1,IPTS: The points at which solution values are
C                         required.
C       UP(K,J)         : Empty array that holds the solution values.
C              K = 1,NPDE AND J = 1,IPTS
C       Outputs
C       *******
C       C(I,J)           J = 1,3 I = 1, N-1 = the polynomial coeffs of
C                        the cubic spline interpolation with interior
C                        knots X(2) to X(N-1). In the interval
C                        X(I)..X(I+1) the spline is given by
C                F(X) = U(K,I) + H* (C(I,1) + H*(C(I,2) +H*C(I,4)))
C                        where H = X -X(I).
C                N.B.    IN THE CASE WHEN NPDE > 1 THESE COEFFICIENTS
C                        WILL BE THOSE FOR TH NPDE TH P.D.E.
C
C       UP(K,JJ)         K = 1,NPDE . JJ = 1,NP .
C                        Array holding solution values for K PDEs at
C                        the points defined by XP(JJ) , JJ =1,NP.
C       IFAIL            Error indicator that is set to 1 if
C                        any of the points XP(JJ) lie outside the
C                        range X(1),X(NPTS).
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IPTS, NPDE, NPTS
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NPTS,3), U(NPDE,NPTS), UP(NPDE,IPTS), X(NPTS),
     *                  XP(IPTS)
C     .. Scalars in Common ..
      DOUBLE PRECISION  DUNFLO, UROUND
      INTEGER           IOVFLO
C     .. Local Scalars ..
      DOUBLE PRECISION  DIVDF1, DIVDF3, DX, P
      INTEGER           I, J, JI, K, KK, L, M, N
      CHARACTER*200     ERRMSG
C     .. External Subroutines ..
      EXTERNAL          D02NNQ
C     .. Common blocks ..
      COMMON            /CD03PC/DUNFLO, UROUND, IOVFLO
C     .. Save statement ..
      SAVE              /CD03PC/
C     .. Executable Statements ..
      IF (NPTS.LT.2) THEN
         IFAIL = 1
         ERRMSG =
     *'  D03PRT- cubic spline routine does not work with
     *   NPTS less than 2. '
         CALL D02NNQ(ERRMSG,1,0,0,0,0,0.0D0,0.0D0)
         RETURN
      END IF
      N = NPTS
      L = N - 1
      DO 180 K = 1, NPDE
C        Form spline coeffs for Kth PDE
C        -------------------------------------------------------------
C        | A tridiagonal linear system for the slopes S(I) of U at X(I)
C        | I = 1,N is generated and then solved by gaussian elimination
C        | with S(I) ending up at C(2,I) for all I = 1,N.
C        | The precise form of the system is
C        |  First equation  not a knot condition at l.h.s.
C        |  --------------
C        |    DX(3) * S(1) + (DX(2) + DX(3)) * S(2) = R(1)
C        | where
C        |    R(1) = (DX(3)*DY(2)*(3DX(2)+2DX(3)) + DY(3)*DX(2)**2 ) /
C        |            (DX(2) + DX(3))
C        |
C        |  Middle equations
C        |  ----------------
C        | DX(I+1)* S(I-1) + (DX(I+1)+DX(I))* S(I) + DX(I)* S(I+1)=R(I)
C        | WHERE
C        |    R(I) = 3*(DX(I+1)*DY(I) + DX(I)*DY(I+1)) , I= 2,...,N-1
C        |
C        |  Last equation  not a knot condition at r.h.s.
C        |  -------------
C        | (DX(N)+DX(N-1) * S(N-1) + DX(N-1) * S(N) = R(N)
C        |  where
C        |    R(N) = (DX(N-1) * DY(N) * (3*DX(N) + 2*DX(N-1)) +
C        |                      DY(N-1) * DX(N) * DX(N)) /(DX(N)+DX(N-1)
C        |  and
C        |    DX(I) = X(I)-X(I-1) , DY(I) = (C(I,1)-C(I-1,1))/DX(I)
C        |                I = 2, ... , N.
C        |
C        |  The upper diagonal of the matrix is stored in C(I,3) and
C        |  the righthand side of the system is stored in C(I,2) and is
C        |  overwritten by the solution S(I) that is the deriv of the
C        |  spline interpolant at the point X(I). The arrays C(I,3)
C        |  and C(.,4) are used initially for tempoary storage of the
C        |  differences of the mesh and of first divided differences of
C        |  the solution.
C        |--------------------------------------------------------------
         DO 20 M = 2, NPTS
            C(M,2) = X(M) - X(M-1)
            C(M,3) = (U(K,M)-U(K,M-1))/C(M,2)
   20    CONTINUE
C
C        Form the top row of the matrix
         C(1,1) = (C(3,2)*C(2,3)*(3.0D0*C(2,2)+2.0D0*C(3,2))+C(3,3)
     *            *C(2,2)**2)/((C(2,2)+C(3,2))*C(3,2))
         C(1,2) = C(2,2)/C(3,2) + 1.0D0
C
         DO 40 I = 2, L
C           Perform the inner elimination loop
            P = 2.0D0*(C(I,2)+C(I+1,2)) - C(I-1,2)*C(I+1,2)
            C(I,1) = (3.0D0*(C(I+1,2)*C(I,3)+C(I,2)*C(I+1,3))-C(I-1,1)
     *               *C(I+1,2))/P
            C(I,2) = C(I,2)/P
C           C(I,2) No longer contains X(I) - X(I-1)
   40    CONTINUE
C
C        Last equation for R hand not a knot condition.
         DX = X(N-1) - X(N-2)
         P = DX - C(N-1,2)*(X(N)-X(N-2))
         C(N,1) = ((DX*C(N,3)*(3.0D0*C(N,2)+2.0D0*DX)+C(N-1,3)*C(N,2)
     *            **2)/(C(N,2)+DX)-C(N-1,1)*(C(N,2)+DX))/P
C         --------------------------------------------------------
C         | The equations now read                               |
C         |       S(I) + C(I,2)*S(I+1) = C(I,1) , I = 1,...,N-1  |
C         | and                 S(N)   = C(N,1) .                |
C         | so carry out back substitution                       |
C         --------------------------------------------------------
         DO 60 I = L, 1, -1
            C(I,1) = C(I,1) - C(I,2)*C(I+1,1)
   60    CONTINUE
C         -------------------------------------------------------------
C         | Generate cubic coeffs in each interval i.e. the derivs.
C         -------------------------------------------------------------
         DO 80 I = 2, NPTS
            DX = X(I) - X(I-1)
            DIVDF1 = (U(K,I)-U(K,I-1))/DX
            DIVDF3 = C(I-1,1) + C(I,1) - 2.0*DIVDF1
            C(I-1,2) = (DIVDF1-C(I-1,1)-DIVDF3)/DX
            C(I-1,3) = DIVDF3/DX**2
   80    CONTINUE
         DO 160 JI = 1, IPTS
C           -----------------------------------------------
C           |Search for the grid points that bound XP(JI) |
C           -----------------------------------------------
            DX = XP(JI)
            IF (DX.LT.(X(1)-UROUND) .OR. DX.GT.(X(NPTS)+UROUND)) THEN
C              DX IS OUTSIDE THE MESH
               IFAIL = 1
            END IF
            DO 100 I = 2, NPTS
               J = I - 1
               IF (DX.LT.X(I)) GO TO 120
  100       CONTINUE
            J = NPTS
C           -----------------------------------------------
C           | XP(JI) lies in the interval (X(J) , X(J+1)) |
C           | Set UP(K,JI) using the spline coeffs in R.  |
C           -----------------------------------------------
  120       DX = DX - X(J)
            UP(K,JI) = C(J,3)
            DO 140 KK = 2, 1, -1
               UP(K,JI) = UP(K,JI)*DX + C(J,KK)
  140       CONTINUE
            UP(K,JI) = UP(K,JI)*DX + U(K,J)
  160    CONTINUE
  180 CONTINUE
      RETURN
      END
