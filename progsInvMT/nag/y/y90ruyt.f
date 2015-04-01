      SUBROUTINE Y90RUY(N,NNZ,P,Q,C,S,A,IROW,ICOL,IDIMA,ISTR,ISTC,AP,AQ)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
*-----------------------------------------------------------------------
*
*         ====================================
*         *  Y90RUY :  Auxiliary for Y90RUF  *
*         ====================================
*
*     Purpose
*     =======
*     Calculates an elementary plane rotation on a real symmetric matrix
*     stored in linked-list format.
*
*-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  C, S
      INTEGER           IDIMA, N, NNZ, P, Q
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IDIMA), AP(N), AQ(N)
      INTEGER           ICOL(IDIMA), IROW(IDIMA), ISTC(N), ISTR(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  APP, APQ, AQP, AQQ
C     .. External Subroutines ..
      EXTERNAL          F06EPF, Y90RUW, Y90RUX
C     .. Executable Statements ..
*
*     Copy non-zero elements of row P into AP and row Q into AQ.
*
      CALL Y90RUX(N,P,A,IROW,ICOL,IDIMA,ISTR,ISTC,AP)
      CALL Y90RUX(N,Q,A,IROW,ICOL,IDIMA,ISTR,ISTC,AQ)
*
*     Compute plane rotation of rows P and Q.
*
      CALL F06EPF(N,AP,1,AQ,1,C,S)
*
      APP = C*AP(P) + S*AP(Q)
      APQ = C*AP(Q) - S*AP(P)
      AP(P) = APP
      AP(Q) = APQ
      AQP = C*AQ(P) + S*AQ(Q)
      AQQ = C*AQ(Q) - S*AQ(P)
      AQ(P) = AQP
      AQ(Q) = AQQ
*
*     Map back to linked-list format adding fill elements.
*
      AP(Q) = ZERO
      CALL Y90RUW(N,P,NNZ,A,IROW,ICOL,IDIMA,ISTR,ISTC,AP)
      CALL Y90RUW(N,Q,NNZ,A,IROW,ICOL,IDIMA,ISTR,ISTC,AQ)
*
*     End of subroutine Y90RUX
*
      RETURN
      END
