      SUBROUTINE D01AHV(A,B,AMAXL,AMAXR)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     CALCULATES THE RELATIVE GRADIENTS AT A AND B (RESPECTIVE SIZES
C     AMAXL AND AMAXR)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, AMAXL, AMAXR, B
C     .. Scalars in Common ..
      DOUBLE PRECISION  FZERO
C     .. Arrays in Common ..
      DOUBLE PRECISION  FUNCTM(127), FUNCTP(127)
C     .. Local Scalars ..
      DOUBLE PRECISION  D1, D30, P, Q, SUM, T
      INTEGER           I, J
C     .. Local Arrays ..
      DOUBLE PRECISION  PIV(15)
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /BD01AH/FUNCTP, FUNCTM, FZERO
C     .. Data statements ..
C
C     5 DIGITS ARE SUFFICIENT HERE FOR PIV (31-POINT RULE NODES)
      DATA              PIV(1), PIV(2), PIV(3), PIV(4), PIV(5), PIV(6),
     *                  PIV(7), PIV(8), PIV(9), PIV(10), PIV(11),
     *                  PIV(12), PIV(13), PIV(14), PIV(15)/.99910D0,
     *                  .99383D0, .98153D0, .96049D0, .92965D0,
     *                  .88846D0, .83673D0, .77460D0, .70250D0,
     *                  .62110D0, .53132D0, .43424D0, .33114D0,
     *                  .22339D0, .11249D0/
C     .. Executable Statements ..
      SUM = 0.0D0
      DO 20 J = 1, 14
         I = J*8
         T = PIV(J) - PIV(J+1)
         P = (FUNCTM(I+8)-FUNCTM(I))/T
         Q = (FUNCTP(I)-FUNCTP(I+8))/T
         SUM = SUM + ABS(P) + ABS(Q)
   20 CONTINUE
      T = PIV(15)
      P = (-FUNCTM(120)+FZERO)/T
      Q = (-FZERO+FUNCTP(120))/T
      T = PIV(1) - PIV(2)
      D1 = (FUNCTM(16)-FUNCTM(8))/T
      D30 = (FUNCTP(8)-FUNCTP(16))/T
      SUM = SUM + ABS(P) + ABS(Q)
      AMAXL = ABS(D1)/SUM
      AMAXR = ABS(D30)/SUM
      RETURN
      END
