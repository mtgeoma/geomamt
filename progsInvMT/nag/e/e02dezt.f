      SUBROUTINE E02DEZ(M,NXKNTS,XMIN,XMAX,XLAM,LXLAM,NYKNTS,YMIN,YMAX,
     *                  YLAM,LYLAM,C,IC1,IC2,LC,X,Y,S,YDEP,IYBSI)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 17 REVISED. IER-1635 (JUN 1995).
C     Chopped from DASL routine B2VR0.
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XMAX, XMIN, YMAX, YMIN
      INTEGER           IC1, IC2, LC, LXLAM, LYLAM, M, NXKNTS, NYKNTS
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LC), S(M), X(M), XLAM(LXLAM), Y(M),
     *                  YDEP(LYLAM-4), YLAM(LYLAM)
      INTEGER           IYBSI(LYLAM-4)
C     .. Local Scalars ..
      DOUBLE PRECISION  T
      INTEGER           IC, IND, IXBAS, IYBAS, J, JSTOP, JSTRT, JXINT,
     *                  JYINT, K, NINDEX
C     .. Local Arrays ..
      DOUBLE PRECISION  XBASIS(4), YBASIS(4)
      INTEGER           IYINT(1)
C     .. External Subroutines ..
      EXTERNAL          E01DAV, E02DFV, E02DFX
C     .. Executable Statements ..
      JXINT = 0
      DO 80 K = 1, M
C        Determine (i) the set of interval indices
C        corresponding to the specified Y-value, and (ii) the
C        indices of B-splines which contain within their
C        support the intervals specified by the interval
C        indices in (i), and the number of such indices.
         CALL E02DFX(4,NYKNTS+4,YLAM(5),LYLAM-4,YMAX,1,Y(K),IYINT,
     *               NINDEX,IYBSI)
         JYINT = IYINT(1)
C        Form the non-zero B-spline basis values (in Y) for
C        each specified  Y-value.
         CALL E01DAV(4,YLAM(JYINT+1),8,Y(K),YBASIS)
C        Determine the interval index corresponding to X(K).
         CALL E02DFV(NXKNTS,XLAM(5),LXLAM-4,XMAX,X(K),JXINT)
C        Form the non-zero B-spline basis values (in X) at X(K).
         CALL E01DAV(4,XLAM(JXINT+1),8,X(K),XBASIS)
C        Form in YDEP all Y-dependent parts of the bivariate spline
C        needed in computing the required value of S(X,Y) at X = X(K).
         DO 40 IND = 1, NINDEX
            J = IYBSI(IND)
            IC = 1 + (JXINT-1)*IC1 + (J-1)*IC2
            T = ZERO
            DO 20 IXBAS = 1, 4
               IC = IC + IC1
               T = T + XBASIS(IXBAS)*C(IC)
   20       CONTINUE
            YDEP(J) = T
   40    CONTINUE
C        Finally form the value of S(X(K),Y(K)).
         IYBAS = 0
         JSTRT = JYINT + 1
         JSTOP = JYINT + 4
         T = ZERO
         DO 60 J = JSTRT, JSTOP
            IYBAS = IYBAS + 1
            T = T + YBASIS(IYBAS)*YDEP(J)
   60    CONTINUE
         S(K) = T
   80 CONTINUE
      RETURN
      END
