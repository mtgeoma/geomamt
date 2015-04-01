      SUBROUTINE F04MCF(N,L,LL,D,NROW,P,B,NRB,ISELCT,X,NRX,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-746 (DEC 1989).
C
C     ******************************************************
C
C     NPL ALGORITHMS LIBRARY ROUTINE VBSOLC
C
C     CREATED 16 05 79.  UPDATED 07 09 79.  RELEASE 00/12
C
C     AUTHOR ... MAURICE G. COX.
C     NATIONAL PHYSICAL LABORATORY, TEDDINGTON,
C     MIDDLESEX TW11 OLW, ENGLAND.
C
C     ******************************************************
C
C     F04MCF.  AN ALGORITHM TO DETERMINE THE SOLUTION OF THE
C     SYSTEM  AX = B  (AND RELATED SYSTEMS),  WHERE  A  IS A
C     SYMMETRIC POSITIVE DEFINITE VARIABLE BANDWIDTH MATRIX,
C     FOLLOWING AN  LDU  FACTORIZATION OF  A  (WHERE
C     U = L TRANSPOSED).
C
C     INPUT PARAMETERS
C        N        ORDER OF  L
C        L        ELEMENTS WITHIN ENVELOPE OF UNIT LOWER
C                    TRIANGULAR MATRIX  L  (INCLUDING UNIT
C                    DIAGONAL ELEMENTS), IN ROW BY ROW ORDER
C        LL       DIMENSION OF  L.  AT LEAST
C                    NROW(1) + NROW(2) + ... + NROW(N).
C        D        DIAGONAL ELEMENTS OF DIAGONAL MATRIX  D
C        NROW     WIDTHS OF ROWS OF  L
C        P        NUMBER OF COLUMNS OF  B.
C                    GREATER THAN OR EQUAL TO  1.
C        B        RIGHT HAND SIDE VECTORS, STORED BY COLUMNS
C        NRB      ROW DIMENSION OF  B,  AS DECLARED IN CALLING
C                    PROGRAM UNIT
C        ISELCT   SELECTION INDICATOR
C                    1 - SOLVE  LDUX = B  (NORMAL APPLICATION)
C                    2 - SOLVE  LDX  = B  (LOWER TRIANGULAR
C                                            SOLVER)
C                    3 - SOLVE  DUX  = B  (UPPER TRIANGULAR
C                                            SOLVER)
C                    4 - SOLVE  LUX  = B
C                    5 - SOLVE  LX   = B  (UNIT LOWER
C                                            TRIANGULAR SOLVER)
C                    6 - SOLVE  UX   = B  (UNIT UPPER
C                                            TRIANGULAR SOLVER),
C                              WHERE  U = L TRANSPOSED.
C
C     OUTPUT (AND ASSOCIATED) PARAMETERS
C        X        SOLUTION VECTORS, STORED BY COLUMNS
C        NRX      ROW DIMENSION OF  X,  AS DECLARED IN CALLING
C                    PROGRAM UNIT
C
C     FAILURE INDICATORS
C        IFAIL    FAILURE INDICATOR
C                    0 - SUCCESSFUL TERMINATION
C                    1 - AT LEAST ONE OF THE FOLLOWING
C                        RESTRICTIONS RELATING TO  N,  LL
C                        AND THE ARRAY  NROW  IS VIOLATED.
C                           N  .GE. 1,
C                           NROW(1) .GE. 1  AND  .LE. I,
C                              FOR  I = 1, 2, ..., N,
C                           NROW(1) + NROW(2) + ...
C                           ... + NROW(N)  .LE.  LL.
C                    2 - AT LEAST ONE OF THE FOLLOWING
C                        RESTRICTIONS RELATING TO  P,
C                        NRB  AND  NRX  IS VIOLATED.
C                           P .GE. 1,
C                           NRB .GE. N,
C                           NRX .GE. N
C                    3 - ISELCT IS LESS THAN 1 OR GREATER THAN 6
C                    4 - D  SINGULAR
C                           (APPLIES ONLY IF  ISELCT .LE. 3)
C                    5 - L  NOT  U N I T  TRIANGULAR, AS EXPECTED
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04MCF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, ISELCT, LL, N, NRB, NRX, P
C     .. Array Arguments ..
      DOUBLE PRECISION  B(NRB,*), D(*), L(LL), X(NRX,*)
      INTEGER           NROW(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ONE, ZERO
      INTEGER           I, IERROR, K, NL
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F04MCZ
C     .. Data statements ..
      DATA              ZERO, ONE/0.0D+0, 1.0D+0/
C     .. Executable Statements ..
C
C     DATA VALIDATION
C
C     CHECK THAT THE ORDER OF  L  IS AT LEAST UNITY,
C     THAT EACH VALUE OF  NROW(I)  LIES BETWEEN
C     APPROPRIATE LIMITS, AND THAT THE SUM OF VALUES OF  NROW(I),
C     FOR  I = 1, 2, ..., N,  IS NO GREATER THAN THE
C     DIMENSION OF  L
C
      IERROR = 1
      IF (N.LT.1) GO TO 100
      NL = 0
      DO 20 I = 1, N
         NL = NL + NROW(I)
         IF (NROW(I).LT.1 .OR. NROW(I).GT.I) GO TO 100
   20 CONTINUE
      IF (NL.GT.LL) GO TO 100
C
C     CHECK THAT THERE IS AT LEAST ONE RIGHT HAND SIDE VECTOR,
C     AND THAT THE ROW DIMENSIONS OF THE ARRAYS  B  AND  X,  AS
C     DECLARED IN THE CALLING PROGRAM UNIT, ARE COMPATIBLE WITH
C     THE ACTUAL NUMBER OF ROWS OF  B  AND  X,  RESPECTIVELY
C
      IERROR = 2
      IF (P.LT.1 .OR. NRB.LT.N .OR. NRX.LT.N) GO TO 100
C
C     CHECK THAT THE SELECTION INDICATOR HAS A PERMITTED VALUE
C
      IERROR = 3
      IF (ISELCT.LT.1 .OR. ISELCT.GT.6) GO TO 100
C
C     CHECK THAT  D  IS NON-SINGULAR (ONLY IF  ISELCT .LE. 3)
C
      IF (ISELCT.GT.3) GO TO 60
      IERROR = 4
      DO 40 I = 1, N
         IF (D(I).EQ.ZERO) GO TO 100
   40 CONTINUE
C
C     CHECK THAT  L  IS  U N I T  TRIANGULAR
C
   60 IERROR = 5
      K = 0
      DO 80 I = 1, N
         K = K + NROW(I)
         IF (L(K).NE.ONE) GO TO 100
   80 CONTINUE
C
C     SOLVE THE SPECIFIED SYSTEM
C
      CALL F04MCZ(N,L,LL,D,NROW,P,B,1,NRB,NRB*P,ISELCT,X,1,NRX,NRX*P)
C
C     RECORD SUCCESSFUL TERMINATION
C
      IERROR = 0
C
C     SET FAILURE INDICATOR
C
  100 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
