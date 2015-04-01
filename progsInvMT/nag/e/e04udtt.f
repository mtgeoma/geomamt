      SUBROUTINE E04UDT(INFORM,N,NCLIN,NCNLN,ALFA,ALFMIN,ALFMAX,BIGBND,
     *                  DXNORM,ANORM,ADX,AX,BL,BU,DSLK,DX,SLK,X)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 16 REVISED. IER-1096 (JUL 1993).
C     MARK 17 REVISED. IER-1616 (JUN 1995).
C
C     ******************************************************************
C     E04UDT  finds a step ALFA such that the point X + ALFA*P reaches
C     one of the slacks or linear constraints.  The step ALFA is the
C     maximum step that can be taken without violating one of the slacks
C     or linear constraints that is currently satisfied.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original Fortran 77 version written  June 1986.
C     This version of E04UDT dated  13-Jun-1987.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, ALFMAX, ALFMIN, BIGBND, DXNORM
      INTEGER           INFORM, N, NCLIN, NCNLN
C     .. Array Arguments ..
      DOUBLE PRECISION  ADX(*), ANORM(*), AX(*), BL(*), BU(*), DSLK(*),
     *                  DX(N), SLK(*), X(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPSPT3, EPSPT5, EPSPT8, EPSPT9
C     .. Local Scalars ..
      DOUBLE PRECISION  ADXI, AXI, RES, ROWNRM
      INTEGER           I, J
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Common blocks ..
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
C     .. Executable Statements ..
C
      ALFA = ALFMAX
      J = 1
C
C     +    WHILE (J .LE. N+NCLIN+NCNLN .AND. ALFA .GT. ALFMIN) DO
   20 IF (J.LE.N+NCLIN+NCNLN .AND. ALFA.GT.ALFMIN) THEN
C
         IF (J.LE.N) THEN
            AXI = X(J)
            ADXI = DX(J)
            ROWNRM = ONE
         ELSE IF (J.LE.N+NCLIN) THEN
C
            I = J - N
            AXI = AX(I)
            ADXI = ADX(I)
            ROWNRM = ANORM(I) + ONE
         ELSE
C
            I = J - N - NCLIN
            AXI = SLK(I)
            ADXI = DSLK(I)
            ROWNRM = ONE
         END IF
C
         RES = -ONE
         IF (ADXI.LE.-EPSPT9*ROWNRM*DXNORM) THEN
C
C           Constraint decreasing.
C
            ADXI = -ADXI
            IF (BL(J).GT.-BIGBND) RES = AXI - BL(J)
         ELSE IF (ADXI.GT.EPSPT9*ROWNRM*DXNORM) THEN
C
C           Constraint increasing.
C
            IF (BU(J).LT.BIGBND) RES = BU(J) - AXI
         END IF
C
         IF (RES.GT.ZERO .AND. ALFA*ADXI.GT.RES) ALFA = RES/ADXI
C
         J = J + 1
         GO TO 20
C        +    END WHILE
      END IF
C
C     ==================================================================
C     Determine ALFA, the bound on the step to be taken.
C     ==================================================================
      ALFA = MAX(ALFA,ALFMIN)
C
      INFORM = 0
      IF (ALFA.GE.ALFMAX) INFORM = 1
C
      RETURN
C
C
C     End of  E04UDT. (NPALF)
C
      END
