      SUBROUTINE F04MCW(N,L,LL,NROW,P,B,IB1,IB2,LB,X,IX1,IX2,LX)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     ******************************************************
C
C     NPL ALGORITHMS LIBRARY ROUTINE VBUSOL
C
C     CREATED 26 06 79.  UPDATED 23 07 79.  RELEASE 00/04
C
C     AUTHOR ... MAURICE G. COX.
C     NATIONAL PHYSICAL LABORATORY, TEDDINGTON,
C     MIDDLESEX TW11 OLW, ENGLAND.
C
C     ******************************************************
C
C     F04MCW.  AN ALGORITHM TO DETERMINE THE SOLUTION OF THE
C     SYSTEM  UX = B,  WHERE  U  IS THE TRANSPOSE OF
C     A UNIT LOWER UNIT TRIANGULAR VARIABLE BANDWIDTH MATRIX.
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
      DOUBLE PRECISION  ONE, XNOW
      INTEGER           I, IL, ILREF, IP, IREV, IX, IXNOW, IXREF, JMAX,
     *                  JREV
C     .. External Subroutines ..
      EXTERNAL          F04MCV
C     .. Data statements ..
      DATA              ONE/1.0D+0/
C     .. Executable Statements ..
C
C     COPY  B  TO  X
C
      CALL F04MCV(N,P,B,IB1,IB2,LB,ONE,X,IX1,IX2,LX)
C
C     DETERMINE NEXT AVAILABLE LOCATION IN  L
C
      ILREF = 1
      DO 20 I = 1, N
         ILREF = ILREF + NROW(I)
   20 CONTINUE
C
C     COUNT OVER RIGHT HAND SIDES
C
      IXREF = 1 + N*IX1 - IX2
      DO 80 IP = 1, P
         IXREF = IXREF + IX2
C
C        SOLVE  UX = B  (U = L TRANSPOSED) FOR
C        CURRENT RIGHT HAND SIDE
C
         IL = ILREF
         I = N + 1
         IXNOW = IXREF
         DO 60 IREV = 1, N
            IXNOW = IXNOW - IX1
            I = I - 1
            XNOW = X(IXNOW)
            IL = IL - 1
            JMAX = NROW(I) - 1
            IF (JMAX.EQ.0) GO TO 60
            IX = IXNOW
            DO 40 JREV = 1, JMAX
               IL = IL - 1
               IX = IX - IX1
               X(IX) = X(IX) - L(IL)*XNOW
   40       CONTINUE
   60    CONTINUE
   80 CONTINUE
      RETURN
      END
