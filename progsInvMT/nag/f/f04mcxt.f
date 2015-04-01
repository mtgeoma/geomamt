      SUBROUTINE F04MCX(N,L,LL,NROW,P,B,IB1,IB2,LB,X,IX1,IX2,LX)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     ******************************************************
C
C     NPL ALGORITHMS LIBRARY ROUTINE VBLSOL
C
C     CREATED 25 06 79.  UPDATED 23 07 79.  RELEASE 00/05
C
C     AUTHOR ... MAURICE G. COX.
C     NATIONAL PHYSICAL LABORATORY, TEDDINGTON,
C     MIDDLESEX TW11 OLW, ENGLAND.
C
C     ******************************************************
C
C     F04MCX.  AN ALGORITHM TO DETERMINE THE SOLUTION OF THE
C     SYSTEM  LX = B,  WHERE  L  IS A UNIT LOWER UNIT TRIANGULAR
C     VARIABLE BANDWIDTH MATRIX.
C     NO CHECKS ARE PERFORMED.
C     B  AND  X  ARE REPRESENTED BY RECTANGULAR DATA STRUCTURES.
C
C     INPUT PARAMETERS
C        N        ORDER OF  L
C        L        ELEMENTS WITHIN ENVELOPE OF UNIT LOWER
C                    TRIANGULAR MATRIX  L  (INCLUDING UNIT
C                    DIAGONAL ELEMENTS), IN ROW BY ROW ORDER
C        LL       DIMENSION OF  L.  AT LEAST
C                    NROW(1) + NROW(2) + ... + NROW(N).
C        NROW     WIDTHS OF ROWS OF  L
C        P        NUMBER OF COLUMNS OF  B.
C                    GREATER THAN OR EQUAL TO  1.
C        B        RIGHT HAND SIDE VECTORS
C        IB1      FIRST INDEX INCREMENT FOR  B
C        IB2      SECOND INDEX INCREMENT FOR  B
C        LB       DIMENSION OF  B
C
C     OUTPUT (AND ASSOCIATED) PARAMETERS
C        X        SOLUTION VECTORS.  MAY BE ASSOCIATED WITH  B.
C        IX1      FIRST INDEX INCREMENT FOR  X
C        IX2      SECOND INDEX INCREMENT FOR  X
C        LX       DIMENSION OF  X
C
C     .. Scalar Arguments ..
      INTEGER           IB1, IB2, IX1, IX2, LB, LL, LX, N, P
C     .. Array Arguments ..
      DOUBLE PRECISION  B(LB), L(LL), X(LX)
      INTEGER           NROW(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  T, ZERO
      INTEGER           I, IB, IBREF, IL, IP, IX, IXREF, J, JMIN
C     .. Data statements ..
      DATA              ZERO/0.0D+0/
C     .. Executable Statements ..
      IBREF = 1 - IB1 - IB2
      IXREF = 1 - IX1 - IX2
C
C     COUNT OVER RIGHT HAND SIDES
C
      DO 80 IP = 1, P
         IBREF = IBREF + IB2
         IXREF = IXREF + IX2
C
C        SOLVE  LX = B  FOR CURRENT RIGHT HAND SIDE
C
         IB = IBREF
         IL = -1
         DO 60 I = 1, N
            IL = IL + 1
            T = ZERO
            JMIN = I - NROW(I)
            IX = IXREF + JMIN*IX1
            JMIN = JMIN + 2
            IF (JMIN.GT.I) GO TO 40
            DO 20 J = JMIN, I
               IL = IL + 1
               IX = IX + IX1
               T = T + L(IL)*X(IX)
   20       CONTINUE
   40       IB = IB + IB1
            IX = IX + IX1
            X(IX) = B(IB) - T
   60    CONTINUE
   80 CONTINUE
      RETURN
      END
