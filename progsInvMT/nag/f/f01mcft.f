      SUBROUTINE F01MCF(N,A,LAL,NROW,L,D,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-729 (DEC 1989).
C
C     ******************************************************
C
C     NPL ALGORITHMS LIBRARY ROUTINE VBCHFC
C
C     CREATED 04 05 79.  UPDATED 07 09 79.  RELEASE 00/22
C
C     AUTHOR ... MAURICE G. COX.
C     NATIONAL PHYSICAL LABORATORY, TEDDINGTON,
C     MIDDLESEX TW11 OLW, ENGLAND.
C
C     ******************************************************
C
C     F01MCF.  AN ALGORITHM TO DETERMINE THE CHOLESKY
C     FACTORIZATION  A = LDU  (WHERE  L  IS UNIT TRIANGULAR,
C     D  IS DIAGONAL AND  U = L TRANSPOSED) OF A SYMMETRIC
C     POSITIVE DEFINITE MATRIX  A  HAVING A VARIABLE BANDWIDTH.
C
C     INPUT PARAMETERS
C        N        ORDER OF  A
C        A        ELEMENTS WITHIN ENVELOPE OF  A,
C                    IN ROW BY ROW ORDER
C        LAL      SMALLER OF THE DIMENSIONS OF A AND L. AT LEAST
C                    NROW(1) + NROW(2) + ... + NROW(N).
C        NROW     WIDTHS OF ROWS OF  A
C
C     OUTPUT (AND ASSOCIATED) PARAMETERS
C        L        ELEMENTS WITHIN ENVELOPE OF UNIT LOWER
C                    TRIANGULAR FACTOR  L  (INCLUDING UNIT
C                    DIAGONAL ELEMENTS), IN ROW BY ROW ORDER
C        D        DIAGONAL ELEMENTS OF D
C
C     INPUT-OUTPUT PARAMETER
C        IFAIL    FAILURE INDICATOR
C
C                 ENTRY VALUES
C                    0 - TERMINATE IMMEDIATELY IF MATRIX FOUND TO BE
C                        NON-DEFINITE
C                    NON-ZERO - CONTINUE UNLESS ZERO PIVOT FOUND
C
C                  EXIT VALUES
C                    0 - SUCCESSFUL TERMINATION
C                    1 - AT LEAST ONE OF THE FOLLOWING
C                        RESTRICTIONS RELATING TO  N,  LAL
C                        AND THE ARRAY  NROW  IS VIOLATED.
C                           N  .GE. 1,
C                           NROW(1) .GE. 1  AND  .LE. I,
C                              FOR  I = 1, 2, ..., N,
C                           NROW(1) + NROW(2) + ...
C                    2 - A, AS MODIFIED BY ROUNDING ERRORS,
C                        IS NOT POSITIVE DEFINITE. FACTORIZATION
C                        NOT COMPLETED.
C                    3 - A, AS MODIFIED BY ROUNDING ERRORS,
C                        IS NOT POSITIVE DEFINITE. FACTORIZATION
C                        COMPLETED BUT MAY BE INACCURATE.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01MCF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LAL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LAL), D(*), L(LAL)
      INTEGER           NROW(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ONE, X, Y, Z, ZERO
      INTEGER           I, ICOL, IERROR, IL, ILREF, J, JCOL, JL, K,
     *                  KMAX, KMIN, NA, PROWI, PROWI1, PROWJ, PROWJ1,
     *                  PTR
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MOD
C     .. Data statements ..
      DATA              ZERO, ONE/0.0D+0, 1.0D+0/
C     .. Executable Statements ..
C
C
C     DATA VALIDATION
C
C     CHECK THAT THE ORDER OF  A  IS AT LEAST UNITY,
C     THAT EACH VALUE OF  NROW(I)  LIES BETWEEN
C     APPROPRIATE LIMITS, AND THAT THE SUM OF VALUES OF  NROW(I),
C     FOR  I = 1, 2, ..., N,  IS NO GREATER THAN THE SMALLER OF THE
C     DIMENSIONS OF  A  AND  L
C
      IERROR = 1
      IF (N.LT.1) GO TO 300
      NA = 0
      DO 20 I = 1, N
         NA = NA + NROW(I)
         IF (NROW(I).LT.1 .OR. NROW(I).GT.I) GO TO 300
   20 CONTINUE
      IF (NA.GT.LAL) GO TO 300
C
C     ASSUME POSITIVE DEFINITENESS UNTIL ESTABLISHED OTHERWISE
C
      IERROR = 0
C
C     OVERWRITE THE VALUES OF  NROW(I)  BY THOSE OF  PROW(I)
C
      NROW(1) = NROW(1) + 1
      IF (N.EQ.1) GO TO 60
      DO 40 I = 2, N
         NROW(I) = NROW(I) + NROW(I-1)
   40 CONTINUE
C
C     PTR  POINTS TO THE ELEMENT IN  L  CURRENTLY BEING FORMED
C
   60 PTR = 0
C
C     PROWI  HOLDS THE CURRENT AND  PROWI1  THE PREVIOUS
C     VALUE OF  PROW(I)
C
      PROWI = 1
      DO 220 I = 1, N
C
C        DETERMINE  I-TH  ROW OF  L  AND  D
C
         PROWI1 = PROWI
         PROWI = NROW(I)
C
C        ILREF  IS REQUIRED SEVERAL TIMES SUBSEQUENTLY
C
         ILREF = PROWI - 2 - I
C
C        ICOL  IS THE COLUMN POSITION OF THE FIRST NON-ZERO
C        ELEMENT IN THE  I-TH  ROW
C
         ICOL = I - PROWI + PROWI1 + 1
C
C        PROWJ  HOLDS THE CURRENT AND  PROWJ1  THE PREVIOUS
C        VALUE OF  PROW(J)
C
         IF (ICOL.EQ.1) PROWJ = 1
         IF (ICOL.GT.1) PROWJ = NROW(ICOL-1)
         DO 200 J = ICOL, I
            PROWJ1 = PROWJ
            PROWJ = NROW(J)
C
C           JCOL  IS THE COLUMN POSITION OF THE FIRST NON-ZERO
C           ELEMENT IN THE  J-TH  ROW
C
            JCOL = J - PROWJ + PROWJ1 + 1
            PTR = PTR + 1
            X = A(PTR)
            IF (I.NE.J) GO TO 140
C
C           I  EQUAL TO  J
C
C           SET  I-TH  UNIT DIAGONAL ELEMENT OF  L
C
            L(PTR) = ONE
C
C           DETERMINE ELEMENTS IN  I-TH  ROW OF  L  FROM THOSE
C           OF  L*D  AND FORM THE  I-TH  DIAGONAL ELEMENT OF  D
C
            KMAX = J - 1
            IF (JCOL.GT.KMAX) GO TO 100
C
C           IL  POINTS TO THE CURRENT ELEMENT IN ROW  I
C
            IL = ILREF + JCOL
            DO 80 K = JCOL, KMAX
               IL = IL + 1
               Y = L(IL)
               Z = Y/D(K)
               L(IL) = Z
               X = X - Y*Z
   80       CONTINUE
C
C           SET ERROR FLAG IF  A,  MODIFIED BY THE ROUNDING ERRORS,
C           IS NOT POSITIVE DEFINITE. EXIT IF HARD FAILURE
C           SPECIFIED BY USER BEFORE ENTRY OR PIVOT IS ZERO.
C
  100       IF (X.GT.ZERO) GO TO 120
            IERROR = 3
            IF (X.LT.ZERO .AND. MOD(IFAIL,10).NE.0) GO TO 120
            IERROR = 2
            GO TO 240
  120       D(I) = X
            GO TO 200
C
C           I  NOT EQUAL TO  J
C
C           DETERMINE  J-TH  ELEMENT WITHIN  I-TH  ROW OF  L*D
C
  140       KMIN = MAX(ICOL,JCOL)
            KMAX = J - 1
            IF (KMIN.GT.KMAX) GO TO 180
C
C           JL  POINTS TO THE CURRENT ELEMENT IN ROW  J
C
            IL = ILREF + KMIN
            JL = PROWJ - 2 - J + KMIN
            DO 160 K = KMIN, KMAX
               IL = IL + 1
               JL = JL + 1
               X = X - L(IL)*L(JL)
  160       CONTINUE
  180       L(PTR) = X
  200    CONTINUE
  220 CONTINUE
C
C     RESTORE THE SPECIFIED VALUES OF  NROW(I)
C
  240 IF (N.EQ.1) GO TO 280
      I = N + 1
      DO 260 J = 2, N
         I = I - 1
         NROW(I) = NROW(I) - NROW(I-1)
  260 CONTINUE
  280 NROW(1) = NROW(1) - 1
C
C     SET FAILURE INDICATOR
C
  300 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
