      SUBROUTINE D02PDT(V,HAVG,X,Y,F,FXY,WT,SCALE,VDOTV,Z,ZDOTZ,VTEMP)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C$$$$ SUBROUTINE STIFFD $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C************************************************
C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
C************************************************
C
C  External:     F
C  Input:        V(*), HAVG, X, Y(*), FXY(*), WT(*), SCALE, VDOTV,
C  Output:       Z(*), ZDOTZ
C  Workspace:    VTEMP(*)
C
C  For an input vector V(*) of length NEQ, this subroutine computes a
C  vector Z(*) that approximates the product HAVG*J*V where HAVG is an
C  input scalar and J is the Jacobian matrix of a function F evaluated
C  at the input arguments (X,Y(*)). This function is defined by a
C  subroutine of the form F(T,U,F) that when given T and U(*), returns
C  the value of the function in F(*). The input vector FXY(*) is defined
C  by F(X,Y,FXY). Scaling is a delicate matter. A weighted Euclidean
C  norm is used with the (positive) weights provided in WT(*). The input
C  scalar SCALE is the square root of the unit roundoff times the norm
C  of Y(*). The square of the norm of the input vector V(*) is input as
C  VDOTV. The routine outputs the square of the norm of the output
C  vector Z(*) as ZDOTZ. The subroutine calls the DOUBLE PRECISION
C  FUNCTION D02PDQ(U,V,WT,NEQ) to compute the dot (inner) product. The
C  vector VTEMP(*) is used for working storage.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HAVG, SCALE, VDOTV, X, ZDOTZ
C     .. Array Arguments ..
      DOUBLE PRECISION  FXY(*), V(*), VTEMP(*), WT(*), Y(*), Z(*)
C     .. Subroutine Arguments ..
      EXTERNAL          F
C     .. Scalars in Common ..
      DOUBLE PRECISION  DIR, H, HOLD, HSTRT, T, TND, TOLD, TOLR, TSTRT
      INTEGER           FLSTP, NEQN, NFCN, OKSTP, SVNFCN
      LOGICAL           FIRST, LAST
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP1, TEMP2
      INTEGER           L
C     .. External Functions ..
      DOUBLE PRECISION  D02PDQ
      EXTERNAL          D02PDQ
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Common blocks ..
      COMMON            /AD02PD/TSTRT, TND, DIR, HSTRT, TOLR, NEQN
      COMMON            /BD02PD/T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP,
     *                  FLSTP, FIRST, LAST
C     .. Save statement ..
      SAVE              /AD02PD/, /BD02PD/
C     .. Executable Statements ..
C
C  Scale V(*) so that it can be used as an increment to Y(*)
C  for an accurate difference approximation to the Jacobian.
C
      TEMP1 = SCALE/SQRT(VDOTV)
      DO 20 L = 1, NEQN
         VTEMP(L) = Y(L) + TEMP1*V(L)
   20 CONTINUE
C
      CALL F(X,VTEMP,Z)
      NFCN = NFCN + 1
C
C  Form the difference approximation.  At the same time undo
C  the scaling of V(*) and introduce the factor of HAVG.
C
      TEMP2 = HAVG/TEMP1
      DO 40 L = 1, NEQN
         Z(L) = TEMP2*(Z(L)-FXY(L))
   40 CONTINUE
C
      ZDOTZ = D02PDQ(Z,Z,WT,NEQN)
C
      RETURN
      END
