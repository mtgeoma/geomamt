      SUBROUTINE F04MCY(N,D,P,B,IB1,IB2,LB,X,IX1,IX2,LX)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     ******************************************************
C
C     NPL ALGORITHMS LIBRARY ROUTINE DSOL
C
C     CREATED 26 06 79.  UPDATED 23 07 79.  RELEASE 00/04
C
C     AUTHOR ... MAURICE G. COX.
C     NATIONAL PHYSICAL LABORATORY, TEDDINGTON,
C     MIDDLESEX TW11 OLW, ENGLAND.
C
C     ******************************************************
C
C     SOLVE DIAGONAL SYSTEM, WHERE THE RIGHT HAND SIDES AND
C     SOLUTION VECTORS ARE SPECIFIED BY RECTANGULAR DATA
C     STRUCTURES.
C
C     INPUT PARAMETERS
C        N         ORDER OF MATRIX  D
C        D         DIAGONAL ELEMENTS OF  D
C        P         NUMBER OF RIGHT HAND SIDES
C        B         RIGHT HAND SIDES
C        IB1       FIRST INDEX INCREMENT FOR  B
C        IB2       SECOND INDEX INCREMENT FOR  B
C        LB        DIMENSION OF B
C
C     OUTPUT (AND ASSOCIATED) PARAMETERS
C        X         SOLUTION VECTORS.  MAY BE ASSOCIATED WITH  B.
C        IX1       FIRST INDEX INCREMENT FOR  X
C        IX2       SECOND INDEX INCREMENT FOR  X
C        LX        DIMENSION OF  X
C
C
C     .. Scalar Arguments ..
      INTEGER           IB1, IB2, IX1, IX2, LB, LX, N, P
C     .. Array Arguments ..
      DOUBLE PRECISION  B(LB), D(N), X(LX)
C     .. Local Scalars ..
      INTEGER           IB, IBREF, ID, IP, IX, IXREF
C     .. Executable Statements ..
      IBREF = 1 - IB1 - IB2
      IXREF = 1 - IX1 - IX2
C
C     COUNT OVER RIGHT HAND SIDES
C
      DO 40 IP = 1, P
         IBREF = IBREF + IB2
         IXREF = IXREF + IX2
         IB = IBREF
         IX = IXREF
C
C        COUNT WITHIN CURRENT RIGHT HAND SIDE
C
         DO 20 ID = 1, N
            IB = IB + IB1
            IX = IX + IX1
            X(IX) = B(IB)/D(ID)
   20    CONTINUE
   40 CONTINUE
      RETURN
      END
