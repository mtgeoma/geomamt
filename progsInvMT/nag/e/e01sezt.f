      SUBROUTINE E01SEZ(N,X,Y,Z,RNQ,FNODES,MINNQ,C,B)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     This subroutine serves to construct quadratic nodal
C     function coefficients from the input data X,Y,Z, using
C     Nag routine F04JGF.
C
C     Routine created - December 1986
C     Author          - Richard Franke
C                       Naval Postgraduate School Monterey,
C                       California  93940
C                       Adapted for Nag by H.Scullion (Leic Univ.)
C                       and I. Gladwell (Nag Ltd.)
C
C     Input Parameters:
C
C           N   -  The number of data points.
C
C       X,Y,Z   -  The data points, (X(I),Y(I),Z(I),I=1,N).
C
C         RNQ   -  The radius for the nodal functions.
C
C     Output Parameters:
C
C      FNODES   -  Real array of dimension at least (5*N).
C                  This array is used to store the coefficients for
C                  the nodal functions.
C
C       MINNQ   -  The smallest number of neighbouring data points
C                  used to define any nodal function.
C
C     Workspace:
C
C           C   -  Real array of dimension at least (5*N).
C                  Workspace for F04JGF.
C
C           B   -  Real array of dimension at least (N).
C                  Workspace for F04JGF.
C
C     .. Parameters ..
      INTEGER           NC
      PARAMETER         (NC=5)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RNQ
      INTEGER           MINNQ, N
C     .. Array Arguments ..
      DOUBLE PRECISION  B(N), C(N,NC), FNODES(NC,N), X(N), Y(N), Z(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  D, RMD, SIGMA, TOL, XD, YD
      INTEGER           I, IRANK, J, KER, NA, NPTREG
      LOGICAL           SVD
C     .. Local Arrays ..
      DOUBLE PRECISION  WRK(4*NC)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          F04JGF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN, SQRT
C     .. Executable Statements ..
C
      DO 80 I = 1, N
         NPTREG = 0
         DO 20 J = 1, N
            IF (I.NE.J) THEN
               XD = X(J) - X(I)
               YD = Y(J) - Y(I)
               D = SQRT(XD**2+YD**2)
               IF (RNQ.GT.D .AND. D.GT.ZERO) THEN
                  RMD = (RNQ-D)/D
                  NPTREG = NPTREG + 1
                  C(NPTREG,1) = RMD*XD
                  C(NPTREG,2) = RMD*YD
                  C(NPTREG,3) = C(NPTREG,1)*XD
                  C(NPTREG,4) = C(NPTREG,1)*YD
                  C(NPTREG,5) = C(NPTREG,2)*YD
                  B(NPTREG) = (Z(J)-Z(I))*RMD
               END IF
            END IF
   20    CONTINUE
C
         MINNQ = MIN(MINNQ,NPTREG)
         DO 40 J = 1, NC
            FNODES(J,I) = ZERO
   40    CONTINUE
         IF (NPTREG.EQ.1) THEN
            D = C(1,1)**2 + C(1,2)**2
            FNODES(1,I) = (C(1,1)/D)*B(1)
            FNODES(2,I) = (C(1,2)/D)*B(1)
         ELSE IF (NPTREG.GE.2) THEN
C           Solve the linear least-squares problem Cx=b.
            NA = NC
            IF (NPTREG.LT.NC) NA = 2
            TOL = SQRT(X02AJF())
            KER = 0
            CALL F04JGF(NPTREG,NA,C,N,B,TOL,SVD,SIGMA,IRANK,WRK,4*NA,
     *                  KER)
            DO 60 J = 1, NA
               FNODES(J,I) = B(J)
   60       CONTINUE
         END IF
   80 CONTINUE
C
      RETURN
      END
