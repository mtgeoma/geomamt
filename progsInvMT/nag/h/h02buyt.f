      SUBROUTINE H02BUY(INPUT,IPRINT,MPSLST,NROWA,MAXN,MAXDIM,N,NCLIN,
     *                  NCTOTL,NMPROB,NMOBJ,NMRHS,NMRNG,NMBND,IOBJ,
     *                  BLDEF,BUDEF,NAMES,A,BL,BU,X,INTVAR,ISTATE,REC,
     *                  IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     ****************************************************************
C     H02BUY READS DATA IN STANDARD MPS FORMAT, SPECIFYING CERTAIN
C     QUANTITIES IN THE FOLLOWING ORDER :-
C
C     ROWS SECTION.
C     CONSTRAINT TYPES AND NAMES.  THE TYPES ALLOWED ARE  E, G, L, N.
C
C     COLUMNS SECTION.
C     COEFFICIENTS FOR THE COLUMNS OF  A.
C     UNSPECIFIED COEFFICIENTS WILL BE SET TO ZERO.
C
C     RHS SECTION.
C     COEFFICIENTS FOR A RIGHT-HAND-SIDE VECTOR  B (BL and BU)
C     UNSPECIFIED COEFFICIENTS WILL BE SET TO ZERO.
C
C     RANGES SECTION.
C     PERMISSIBLE RANGES ON THE RHS VECTOR  B.
C
C     BOUNDS SECTION.
C     LOWER AND UPPER BOUNDS ON THE VARIABLES  X.
C     BOUND TYPES ALLOWED ARE FR, FX, LO, MI, PL, UP, BV, UI, LI.
C     MORE THAN ONE BOUNDS SET MAY APPEAR.
C
C     THE RANGES AND BOUNDS SECTIONS ARE OPTIONAL.
C
C     INPUT PARAMETERS
C     ----------------
C
C     INPUT   IS THE INPUT CHANNEL NUMBER FOR THE MPS DATA.
C     IPRINT  IS THE OUTPUT CHANNEL NUMBER FOR DIAGNOSTICS.
C     MPSLST  IF .TRUE., THEN DATA RECORDS ARE LISTED.
C     NROWA   IS THE DECLARED ROW DIMENSION OF THE ARRAY  A  BELOW.
C     MAXN    IS AN UPPER LIMIT ON THE DIMENSION OF X.
C     MAXDIM  IS AN UPPER LIMIT ON THE DIMENSION OF THE PROBLEM, NAMELY
C             AN OVER-ESTIMATE OF THE NO. OF VARIABLES PLUS CONSTRAINTS.
C     NMOBJ   IS THE PARTICULAR ROW OF  A  TO BE SINGLED OUT AS THE ROW
C             THAT WILL BE POINTED TO BY THE OUTPUT VARIABLE  IOBJ.
C             IF NMOBJ IS INPUT AS BLANK, THE FIRST ROW OF TYPE  N
C             THAT IS ENCOUNTERED DURING INPUT WILL BE SELECTED.
C             THERE NEED NOT BE ANY SUCH ROW.
C     NMRHS   IS THE PARTICULAR  RHS  SET TO BE SELECTED.
C             IF NMRHS IS BLANK, THE FIRST  RHS  SET WILL BE USED.
C     NMRNG   IS THE PARTICULAR RANGE SET TO BE SELECTED.
C             IF NMRNG IS BLANK, THE FIRST RANGE SET WILL BE USED.
C     NMBND   IS THE PARTICULAR BOUND SET TO BE SELECTED.
C             IF NMBND IS BLANK, THE FIRST BOUND SET WILL BE USED.
C     BLDEF   IS THE DEFAULT LOWER BOUND FOR EACH COMPONENT OF  X.
C     BUDEF   IS THE DEFAULT UPPER BOUND FOR EACH COMPONENT OF  X.
C
C
C     OUTPUT PARAMETERS
C     -----------------
C
C     NMPROB  THE PROBLEM NAME GIVEN ON THE NAME CARD.
C     NMOBJ   THE SELECTED ROW OF TYPE  N.
C     NMRHS   THE SELECTED RHS   SET.
C     NMRNG   THE SELECTED RANGE SET.
C     NMBND   THE SELECTED BOUND SET.
C     IOBJ    THE NUMBER OF THE CONSTRAINT NAMED NMOBJ.
C     N       THE ACTUAL NUMBER OF VARIABLES FOUND IN THE
C             COLUMNS SECTION.
C     NCLIN   THE ACTUAL NUMBER OF LINEAR CONSTRAINTS FOUND IN THE
C             ROWS SECTION.
C     NCTOTL  THE TOTAL NUMBER OF VARIABLES AND CONSTRAINTS.
C             (NCTOTL = N + NCLIN)
C     ISTATE  RECORDS THE STATE OF EACH VARIABLE  X  AS SPECIFIED
C             BY THE INITIAL BOUNDS SET.
C             ISTATE(J) = 0  IF  X(J)  WAS SPECIFIED BY AN FX BOUND.
C             ISTATE(J) = 1  IF  X(J) = BL(J).
C             ISTATE(J) = 2  IF  X(J) = BU(J).
C     NAMES   AN ARRAY CONTAINING THE NAMES OF THE VARIABLES AND
C             THE CONSTRAINTS, IN THAT ORDER.  (A8 FORMAT).
C             NAMES ARE 8 CHARACTERS LONG.
C     A       AN ARRAY OF DECLARED DIMENSIONS  A(NROWA,MAXN).
C             THE FIRST  NCLIN  ROWS AND THE FIRST  N  COLUMNS OF  A
C             CONTAIN THE LINEAR CONSTRAINT MATRIX  A.
C     BL,BU   LOWER AND UPPER BOUNDS FOR THE VARIABLES,
C             THE LINEAR CONSTRAINTS, AND THE NONLINEAR CONSTRAINTS,
C             IN THAT ORDER.
C     X       A SET OF INITIAL VALUES FOR THE N VARIABLES  X.
C
C
C     ****************************************************************
C
C     .. Parameters ..
      CHARACTER*32      SPACES
      PARAMETER         (SPACES='                                ')
      CHARACTER         DOLLAR
      PARAMETER         (DOLLAR='$')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BLDEF, BUDEF
      INTEGER           IFAIL, INPUT, IOBJ, IPRINT, MAXDIM, MAXN, N,
     *                  NCLIN, NCTOTL, NROWA
      LOGICAL           MPSLST
      CHARACTER*8       NMBND, NMOBJ, NMPROB, NMRHS, NMRNG
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NROWA,MAXN), BL(MAXDIM), BU(MAXDIM), X(MAXN)
      INTEGER           INTVAR(MAXN), ISTATE(MAXDIM)
      CHARACTER*8       NAMES(MAXDIM)
      CHARACTER*80      REC(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGBND
      INTEGER           I, IFAILR, IMARK, IMOVE, IROW, ITYPE, J, JMARK,
     *                  K, M, N1, NLIST, NUMINT
      LOGICAL           FNDNAM, GOTINT, GOTNM, LTRUE
      CHARACTER*8       COLNAM
      CHARACTER*80      RREC
C     .. Local Arrays ..
      DOUBLE PRECISION  AD(2)
      CHARACTER*80      TREC(4)
C     .. External Functions ..
      LOGICAL           H02BUX
      EXTERNAL          H02BUX
C     .. External Subroutines ..
      EXTERNAL          H02BUW, X04BAF, X04BAY, X04BBF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
C     Initialise variables
      BIGBND = 1.0D+20
      GOTINT = .FALSE.
      NUMINT = 0
      M = 0
      N = 0
      FNDNAM = .FALSE.
      IOBJ = 0
      NLIST = 0
      GOTNM = NMOBJ .NE. '        '
      IMARK = 1
C     -----------------------------------------------
C     input problem name, row names and types.
C     -----------------------------------------------
      IF (MPSLST) THEN
         WRITE (TREC,FMT=99990)
         CALL X04BAY(IPRINT,3,TREC)
      END IF
C
C     look for NAME card
C
   20 IFAILR = 1
      CALL X04BBF(INPUT,RREC,IFAILR)
C     test for failure on reading
      IF (IFAILR.NE.0) GO TO 500
C
      NLIST = NLIST + 1
C     If user has requested it, print data
      IF (MPSLST) CALL X04BAF(IPRINT,' '//RREC)
C
C     Ignore comments
      IF (RREC(1:1).EQ.'*') GO TO 20
C
      IF (RREC(1:4).NE.'NAME') GO TO 480
      NMPROB = RREC(15:22)
C
C     *****************************************************
C     Now look for a ROWS card
   40 IFAILR = 1
      CALL X04BBF(INPUT,RREC,IFAILR)
      IF (IFAILR.NE.0) GO TO 500
C
      NLIST = NLIST + 1
C     If user has requested it, print data
      IF (MPSLST) CALL X04BAF(IPRINT,' '//RREC)
C
C     Ignore comments
      IF (RREC(1:1).EQ.'*') GO TO 40
C
C     Expect rows card
      IF (RREC(1:4).EQ.'ROWS') GO TO 60
C
C     Check to see if it is COLUMNS card
      IF (RREC(1:7).EQ.'COLUMNS') THEN
C
C        COLUMNS card read before any constraints.
C
         IF (IFAIL.EQ.0) THEN
            WRITE (REC,FMT=99986)
            IFAIL = 4
         END IF
         IF ( .NOT. MPSLST) RETURN
         WRITE (TREC,FMT=99986)
         CALL X04BAY(IPRINT,2,TREC)
         GO TO 120
      END IF
C
C     Illegal card preceding rows card
C
      IF ( .NOT. MPSLST) GO TO 480
      IF (IFAIL.EQ.0) THEN
         IFAIL = 11
         WRITE (REC,FMT=99998) NLIST
      END IF
      WRITE (TREC,FMT=99998) NLIST
      CALL X04BAY(IPRINT,2,TREC)
      IF (RREC(1:6).EQ.'ENDATA') GO TO 440
      GO TO 40
C
C     Now read the row names and check if the inequality is valid
C
   60 IFAILR = 1
      CALL X04BBF(INPUT,RREC,IFAILR)
      IF (IFAILR.NE.0) GO TO 500
C
      NLIST = NLIST + 1
C     If user has requested it, print data
      IF (MPSLST) CALL X04BAF(IPRINT,' '//RREC)
C
C     Ignore comments
      IF (RREC(1:1).EQ.'*') GO TO 60
C
C     check for COLUMNS card
      IF (RREC(1:7).EQ.'COLUMNS') GO TO 120
C
C     Check that the row name is valid
      IF ( .NOT. H02BUX(RREC(5:12))) THEN
         IF (IFAIL.EQ.0) THEN
            IFAIL = 6
            WRITE (REC,FMT=99995) RREC(5:12)
         END IF
         IF ( .NOT. MPSLST) RETURN
         WRITE (TREC,FMT=99995) RREC(5:12)
         CALL X04BAY(IPRINT,2,TREC)
         GO TO 60
      END IF
      IF (RREC(1:4).EQ.' L  ' .OR. RREC(1:4).EQ.'  L ') THEN
         ITYPE = 1
      ELSE IF (RREC(1:4).EQ.' G  ' .OR. RREC(1:4).EQ.'  G ') THEN
         ITYPE = -1
      ELSE IF (RREC(1:4).EQ.' E  ' .OR. RREC(1:4).EQ.'  E ') THEN
         ITYPE = 0
      ELSE IF (RREC(1:4).EQ.' N  ' .OR. RREC(1:4).EQ.'  N ') THEN
         ITYPE = 2
C
C        Check for chosen OBJECTIVE function
         IF (IOBJ.LT.1) THEN
            IF ( .NOT. GOTNM) NMOBJ = RREC(5:12)
            IF (RREC(5:12).EQ.NMOBJ) IOBJ = M + 1
         END IF
      ELSE
C        neither inequality or COLUMNS
         IF (IFAIL.EQ.0) THEN
            IFAIL = 5
            WRITE (REC,FMT=99988) RREC(1:4)
         END IF
         IF ( .NOT. MPSLST) RETURN
         WRITE (TREC,FMT=99988) RREC(1:4)
         CALL X04BAY(IPRINT,2,TREC)
         IF (RREC(1:6).EQ.'ENDATA') GO TO 440
         GO TO 60
      END IF
C
C     Now check for duplicate names of rows
C
      IF (M.GT.0) THEN
         CALL H02BUW(M,NAMES,RREC(5:12),1,M,IMARK,IROW)
         IF (IROW.GT.0) THEN
C
C           row name repeated, ignore.
C
C           If user has requested it, print data
            IF (MPSLST) THEN
               WRITE (TREC,FMT=99987) RREC(5:12)
               CALL X04BAY(IPRINT,2,TREC)
            END IF
            GO TO 60
         END IF
      END IF
C
C     Enter new row name
C
      M = M + 1
      IF (M.GT.NROWA) THEN
C
C        too many rows - the following code counts the number
C        of ROWS data in M before returning.
C
   80    IFAILR = 1
         CALL X04BBF(INPUT,RREC,IFAILR)
         IF (IFAILR.NE.0) GO TO 100
C
C        Ignore comments
         IF (RREC(1:1).EQ.'*') GO TO 80
         IF (RREC(1:7).NE.'COLUMNS') THEN
            M = M + 1
            GO TO 80
         END IF
  100    IF (IFAIL.EQ.0) THEN
            IFAIL = 1
            WRITE (REC,FMT=99977) NROWA, M
         END IF
         IF ( .NOT. MPSLST) RETURN
         WRITE (TREC,FMT=99977) NROWA, M
         CALL X04BAY(IPRINT,2,TREC)
         IF (RREC(1:6).EQ.'ENDATA') GO TO 440
         GO TO 60
      END IF
C
      NAMES(M) = RREC(5:12)
      ISTATE(M) = ITYPE
      GO TO 60
C
C     *****************************************************************
C     Have got COLUMNS data
C
C     Check if ROWS data has been found
C
  120 IF (M.LE.0) THEN
         IF (IFAIL.EQ.0) THEN
            IFAIL = 4
            WRITE (REC,FMT=99986)
         END IF
         IF ( .NOT. MPSLST) RETURN
         WRITE (TREC,FMT=99986)
         CALL X04BAY(IPRINT,2,TREC)
      END IF
C
C     Check whether OBJECTIVE has been found
      IF (IOBJ.EQ.0) THEN
C
C        no objective row found
C
         IF (IFAIL.EQ.0) THEN
            IFAIL = 3
            WRITE (REC,FMT=99978)
         END IF
         IF ( .NOT. MPSLST) RETURN
         WRITE (TREC,FMT=99978)
         CALL X04BAY(IPRINT,2,TREC)
         IF (RREC(1:6).EQ.'ENDATA') GO TO 440
      END IF
C     -----------------------------
C     input columns.
C     -----------------------------
      NCLIN = M
      N = 0
      LTRUE = .TRUE.
C
C     move row names to the end to make room for col names
C
      J = MAXDIM - M
      DO 140 I = M, 1, -1
         NAMES(I+J) = NAMES(I)
  140 CONTINUE
C
C     set  N1  to point to first row name
      N1 = MAXDIM - M + 1
      IMARK = N1
      IMOVE = MAXDIM - M
C
C     read a column card and check for validity
C
  160 IFAILR = 1
      CALL X04BBF(INPUT,RREC,IFAILR)
      IF (IFAILR.NE.0) GO TO 500
C
      NLIST = NLIST + 1
C     If user has requested it, print data
      IF (MPSLST) CALL X04BAF(IPRINT,' '//RREC)
C
C     Ignore comments
      IF (RREC(1:1).EQ.'*') GO TO 160
C     check for $
      IF (RREC(15:15).EQ.DOLLAR) GO TO 160
      IF (RREC(40:40).EQ.DOLLAR) RREC(40:71) = SPACES
C
C     Validate data format
      IF (RREC(1:4).EQ.'RHS ') THEN
         GO TO 240
      ELSE
         IF (RREC(13:14).NE.'  ' .OR. RREC(23:24)
     *       .NE.'  ' .OR. RREC(37:39).NE.'   ' .OR. RREC(48:49)
     *       .NE.'  ' .OR. RREC(62:71).NE.'           ' .OR. RREC(1:4)
     *       .NE.'    ') THEN
            IF ( .NOT. MPSLST) GO TO 480
            IF (IFAIL.EQ.0) THEN
               IFAIL = 11
               WRITE (REC,FMT=99998) NLIST
            END IF
            WRITE (TREC,FMT=99998) NLIST
            CALL X04BAY(IPRINT,2,TREC)
            IF (RREC(1:6).EQ.'ENDATA') GO TO 440
            GO TO 160
         END IF
      END IF
C
      IF ( .NOT. H02BUX(RREC(5:12))) THEN
         IF (IFAIL.EQ.0) THEN
            IFAIL = 7
            WRITE (REC,FMT=99997) RREC(5:12)
         END IF
         IF ( .NOT. MPSLST) RETURN
         WRITE (TREC,FMT=99997) RREC(5:12)
         CALL X04BAY(IPRINT,2,TREC)
         IF (RREC(1:6).EQ.'ENDATA') GO TO 440
         GO TO 160
      END IF
C
C     check for integer variable
C
      IF (RREC(15:22).EQ.'''MARKER''') THEN
         IF (RREC(40:47).EQ.'''INTORG''' .AND. GOTINT) THEN
            IF (IFAIL.EQ.0) THEN
               IFAIL = 15
               WRITE (REC,FMT=99996)
            END IF
            IF ( .NOT. MPSLST) RETURN
            WRITE (TREC,FMT=99996)
            CALL X04BAY(IPRINT,2,TREC)
         ELSE IF (RREC(40:47).EQ.'''INTORG''') THEN
            GOTINT = .TRUE.
         ELSE IF ( .NOT. GOTINT .AND. RREC(40:47).EQ.'''INTEND''') THEN
            IF (IFAIL.EQ.0) THEN
               IFAIL = 15
               WRITE (REC,FMT=99996)
            END IF
            IF ( .NOT. MPSLST) RETURN
            WRITE (TREC,FMT=99996)
            CALL X04BAY(IPRINT,2,TREC)
         ELSE IF (RREC(40:47).EQ.'''INTEND''') THEN
            GOTINT = .FALSE.
         END IF
         IF (RREC(1:6).EQ.'ENDATA') GO TO 440
         GO TO 160
      END IF
C
C     Get the numbers
      READ (RREC,FMT=99994,ERR=480) AD(1), AD(2)
C
C     Branch if it is the start of a new column.
C     Otherwise process two entries on each card.
C
      IF (N.GT.0) LTRUE = RREC(5:12) .NE. NAMES(N)
      IF (LTRUE) THEN
C
C        new column
C
         N = N + 1
         IF (N.GT.MAXN) THEN
C
C           Too many columns - following code counts the number
C           of columns N before exiting.
C
            COLNAM = RREC(5:12)
  180       IFAILR = 1
            CALL X04BBF(INPUT,RREC,IFAILR)
            IF (IFAILR.NE.0) GO TO 200
C
            IF (RREC(15:15).EQ.DOLLAR) GO TO 180
C
C           Ignore comments
            IF (RREC(1:1).EQ.'*') GO TO 180
C
            IF (RREC(1:4).EQ.'    ') THEN
               IF (COLNAM.NE.RREC(5:12)) THEN
                  COLNAM = RREC(5:12)
                  N = N + 1
               END IF
               GO TO 180
            END IF
  200       IF (IFAIL.EQ.0) THEN
               IFAIL = 2
               WRITE (REC,FMT=99976) MAXN, N
            END IF
            IF ( .NOT. MPSLST) RETURN
            WRITE (TREC,FMT=99976) MAXN, N
            CALL X04BAY(IPRINT,2,TREC)
            IF (RREC(1:6).EQ.'ENDATA') GO TO 440
            GO TO 240
         END IF
C
         IF (GOTINT) THEN
            IF (INTVAR(N).LT.1) NUMINT = NUMINT + 1
            INTVAR(N) = 1
         END IF
C
         NAMES(N) = RREC(5:12)
         DO 220 I = 1, M
            A(I,N) = 0.0D0
  220    CONTINUE
C        End of action if new column, N.
      END IF
C
C     check for row names
C
      CALL H02BUW(MAXDIM,NAMES,RREC(15:22),N1,MAXDIM,IMARK,IROW)
      IF (IROW.LT.1) THEN
         IF (IFAIL.EQ.0) THEN
            IFAIL = 12
            WRITE (REC,FMT=99985) RREC(15:22)
         END IF
         IF ( .NOT. MPSLST) RETURN
         WRITE (TREC,FMT=99985) RREC(15:22)
         CALL X04BAY(IPRINT,2,TREC)
         GO TO 160
      END IF
C
      IROW = IROW - IMOVE
      A(IROW,N) = AD(1)
C
C     Now second half of card
      IF (RREC(40:47).NE.'        ') THEN
         CALL H02BUW(MAXDIM,NAMES,RREC(40:47),N1,MAXDIM,IMARK,IROW)
         IF (IROW.LT.1) THEN
            IF (IFAIL.EQ.0) THEN
               IFAIL = 12
               WRITE (REC,FMT=99985) RREC(40:47)
            END IF
            IF ( .NOT. MPSLST) RETURN
            WRITE (TREC,FMT=99985) RREC(40:47)
            CALL X04BAY(IPRINT,2,TREC)
            GO TO 160
         END IF
C
         IROW = IROW - IMOVE
         A(IROW,N) = AD(2)
C
      END IF
      GO TO 160
C
C     *****************************************************************
C     Must be the RHS card
C
  240 IF (RREC(1:3).NE.'RHS') THEN
         IF ( .NOT. MPSLST) GO TO 480
         IF (IFAIL.EQ.0) THEN
            IFAIL = 11
            WRITE (REC,FMT=99998) NLIST
         END IF
         WRITE (TREC,FMT=99998) NLIST
         CALL X04BAY(IPRINT,2,TREC)
         IF (RREC(1:6).EQ.'ENDATA') GO TO 440
         GO TO 280
      END IF
C
C     Were there any columns ?
C
      IF (N.LT.1) THEN
         IF (IFAIL.EQ.0) THEN
            IFAIL = 13
            WRITE (REC,FMT=99984)
         END IF
         IF ( .NOT. MPSLST) RETURN
         WRITE (TREC,FMT=99984)
         CALL X04BAY(IPRINT,2,TREC)
         GO TO 280
      END IF
C     ------------------------------
C     input rhs
C     ------------------------------
      NCTOTL = N + NCLIN
C
C     Move row names again to be contiguous with col names
C
      J = N1 - N - 1
      DO 260 I = N1, MAXDIM
         NAMES(I-J) = NAMES(I)
         BL(I-J) = 0.0D0
  260 CONTINUE
C
      N1 = N + 1
      IMARK = N1
      GOTNM = NMRHS .NE. '        '
      FNDNAM = .FALSE.
C
C     read card and see if it is the rhs we want
C
  280 IFAILR = 1
      CALL X04BBF(INPUT,RREC,IFAILR)
      IF (IFAILR.NE.0) GO TO 500
C
      NLIST = NLIST + 1
C     If user has requested it, print data
      IF (MPSLST) CALL X04BAF(IPRINT,' '//RREC)
C
C     Ignore comments
      IF (RREC(1:1).EQ.'*') GO TO 280
C     check for $
      IF (RREC(15:15).EQ.DOLLAR) GO TO 280
      IF (RREC(40:40).EQ.DOLLAR) RREC(40:71) = SPACES
C
      IF (RREC(1:2).EQ.'EN' .OR. RREC(1:2).EQ.'BO' .OR. RREC(1:2)
     *    .EQ.'RA') GO TO 320
      IF (RREC(1:4).NE.'    ') GO TO 300
C
      IF (RREC(13:14).NE.'  ' .OR. RREC(23:24).NE.'  ' .OR. RREC(37:39)
     *    .NE.'   ' .OR. RREC(48:49).NE.'  ' .OR. RREC(62:71)
     *    .NE.'           ') THEN
         IF ( .NOT. MPSLST) GO TO 480
         IF (IFAIL.EQ.0) THEN
            IFAIL = 11
            WRITE (REC,FMT=99998) NLIST
         END IF
         WRITE (TREC,FMT=99998) NLIST
         CALL X04BAY(IPRINT,2,TREC)
         GO TO 280
      END IF
C
      IF ( .NOT. GOTNM) THEN
         NMRHS = RREC(5:12)
         GOTNM = .TRUE.
      END IF
C
      IF (NMRHS.NE.RREC(5:12)) GO TO 280
C
      FNDNAM = .TRUE.
C     Convert numbers
      READ (RREC,FMT=99994,ERR=480) AD(1), AD(2)
C
C     look at first half of record
C
      CALL H02BUW(NCTOTL,NAMES,RREC(15:22),N1,NCTOTL,IMARK,IROW)
      IF (IROW.LT.1) THEN
         IF (IFAIL.EQ.0) THEN
            IFAIL = 12
            WRITE (REC,FMT=99985) RREC(15:22)
         END IF
         IF ( .NOT. MPSLST) RETURN
         WRITE (TREC,FMT=99985) RREC(15:22)
         CALL X04BAY(IPRINT,2,TREC)
         GO TO 280
      END IF
C
C     another rhs element made it
C
      BL(IROW) = AD(1)
C
C     Now examine second half of record
      IF (RREC(40:47).NE.'        ') THEN
C
         CALL H02BUW(NCTOTL,NAMES,RREC(40:47),N1,NCTOTL,IMARK,IROW)
         IF (IROW.LT.1) THEN
            IF (IFAIL.EQ.0) THEN
               IFAIL = 12
               WRITE (REC,FMT=99985) RREC(40:47)
            END IF
            IF ( .NOT. MPSLST) RETURN
            WRITE (TREC,FMT=99985) RREC(40:47)
            CALL X04BAY(IPRINT,2,TREC)
            GO TO 280
         END IF
C
C        another rhs element made it
C
         BL(IROW) = AD(2)
C
      END IF
      GO TO 280
C
C     *****************************************************************
C     comment, RANGES, BOUNDS, ENDATA or garbage
C
  300 IF (RREC(1:2).NE.'BO' .AND. RREC(1:2).NE.'EN' .AND. RREC(1:2)
     *    .NE.'RA') THEN
         IF ( .NOT. MPSLST) GO TO 480
         IF (IFAIL.EQ.0) THEN
            IFAIL = 11
            WRITE (REC,FMT=99998) NLIST
         END IF
         WRITE (TREC,FMT=99998) NLIST
         CALL X04BAY(IPRINT,2,TREC)
         GO TO 280
      END IF
C
C     check for RHSNAM
C
  320 IF ( .NOT. FNDNAM) THEN
         IF (IFAIL.EQ.0) THEN
            IFAIL = 14
            WRITE (REC,FMT=99993) NMRHS
         END IF
         IF ( .NOT. MPSLST) RETURN
         WRITE (TREC,FMT=99993) NMRHS
         CALL X04BAY(IPRINT,2,TREC)
      END IF
C
C     --------------------------
C     input ranges
C     --------------------------
C     set default bounds on constraints.
C
      DO 340 I = 1, M
         J = N + I
         BU(J) = BL(J)
         K = ISTATE(I)
         IF (K.LT.0) BU(J) = BIGBND
         IF (K.GE.1) BL(J) = -BIGBND
         IF (K.EQ.2) BU(J) = BIGBND
  340 CONTINUE
C
      FNDNAM = .TRUE.
      IF (RREC(1:6).NE.'RANGES') GO TO 380
      GOTNM = NMRNG .NE. '        '
      FNDNAM = .FALSE.
C
C     Read card and see if it is the required range record
C
  360 IFAILR = 1
      CALL X04BBF(INPUT,RREC,IFAILR)
      IF (IFAILR.NE.0) GO TO 500
C
      NLIST = NLIST + 1
C     If user has requested it, print data
      IF (MPSLST) CALL X04BAF(IPRINT,' '//RREC)
C
C     Ignore comments
      IF (RREC(1:1).EQ.'*') GO TO 360
C     check for $
      IF (RREC(15:15).EQ.DOLLAR) GO TO 360
      IF (RREC(40:40).EQ.DOLLAR) RREC(40:71) = SPACES
C
      IF (RREC(1:4).NE.'    ') THEN
C        BOUNDS, ENDATA or garbage
         IF (RREC(1:6).NE.'BOUNDS' .AND. RREC(1:6).NE.'ENDATA') THEN
            IF ( .NOT. MPSLST) GO TO 480
            IF (IFAIL.EQ.0) THEN
               IFAIL = 11
               WRITE (REC,FMT=99998) NLIST
            END IF
            WRITE (TREC,FMT=99998) NLIST
            CALL X04BAY(IPRINT,2,TREC)
            GO TO 360
         END IF
C
         GO TO 380
      END IF
C
      IF (RREC(13:14).NE.'  ' .OR. RREC(23:24).NE.'  ' .OR. RREC(37:39)
     *    .NE.'   ' .OR. RREC(48:49).NE.'  ' .OR. RREC(62:71)
     *    .NE.'           ') THEN
         IF ( .NOT. MPSLST) GO TO 480
         IF (IFAIL.EQ.0) THEN
            IFAIL = 11
            WRITE (REC,FMT=99998) NLIST
         END IF
         WRITE (TREC,FMT=99998) NLIST
         CALL X04BAY(IPRINT,2,TREC)
         GO TO 360
      END IF
C
      IF ( .NOT. GOTNM) THEN
         NMRNG = RREC(5:12)
         GOTNM = .TRUE.
      END IF
C
      IF (NMRNG.NE.RREC(5:12)) GO TO 360
C
      FNDNAM = .TRUE.
C     Look at first half of record
C
      CALL H02BUW(NCTOTL,NAMES,RREC(15:22),N1,NCTOTL,IMARK,IROW)
      IF (IROW.LT.1) THEN
         IF (IFAIL.EQ.0) THEN
            IFAIL = 12
            WRITE (REC,FMT=99985) RREC(15:22)
         END IF
         IF ( .NOT. MPSLST) RETURN
         WRITE (TREC,FMT=99985) RREC(15:22)
         CALL X04BAY(IPRINT,2,TREC)
         GO TO 360
      END IF
C
C     Another range element made it
C
C     Convert numbers
      READ (RREC,FMT=99994,ERR=480) AD(1), AD(2)
      J = IROW
      K = ISTATE(IROW-N)
      IF (K.LT.0) BU(J) = BL(J) + ABS(AD(1))
      IF (K.GT.0) BL(J) = BU(J) - ABS(AD(1))
      IF (K.EQ.0 .AND. AD(1).GT.0.0D0) BU(J) = BL(J) + AD(1)
      IF (K.EQ.0 .AND. AD(1).LT.0.0D0) BL(J) = BU(J) + AD(1)
C
C     Examine second half of record
      IF (RREC(40:47).NE.'        ') THEN
         CALL H02BUW(NCTOTL,NAMES,RREC(40:47),N1,NCTOTL,IMARK,IROW)
         IF (IROW.LT.1) THEN
            IF (IFAIL.EQ.0) THEN
               IFAIL = 12
               WRITE (REC,FMT=99985) RREC(40:47)
            END IF
            IF ( .NOT. MPSLST) RETURN
            WRITE (TREC,FMT=99985) RREC(40:47)
            CALL X04BAY(IPRINT,2,TREC)
            GO TO 360
         END IF
C
C        Another range element made it
C
         J = IROW
         K = ISTATE(IROW-N)
         IF (K.LT.0) BU(J) = BL(J) + ABS(AD(2))
         IF (K.GT.0) BL(J) = BU(J) - ABS(AD(2))
         IF (K.EQ.0 .AND. AD(2).GT.0.0D0) BU(J) = BL(J) + AD(2)
         IF (K.EQ.0 .AND. AD(2).LT.0.0D0) BL(J) = BU(J) + AD(2)
C
      END IF
C
      GO TO 360
C
C     *****************************************************************
C     --------------------
C     input bounds.
C     --------------------
C     set bounds to default values
C
  380 DO 400 J = 1, N
         BL(J) = BLDEF
         BU(J) = BUDEF
  400 CONTINUE
C
C     check for RANGE name
      IF ( .NOT. FNDNAM) THEN
         IF (IFAIL.EQ.0) THEN
            IFAIL = 14
            WRITE (REC,FMT=99992) NMRNG
         END IF
         IF ( .NOT. MPSLST) RETURN
         WRITE (TREC,FMT=99992) NMRNG
         CALL X04BAY(IPRINT,2,TREC)
      END IF
C
C     check for no bounds
C
      IF (RREC(1:6).EQ.'ENDATA') GO TO 440
      JMARK = 1
      GOTNM = NMBND .NE. '        '
      FNDNAM = .FALSE.
C
C     read and check bounds cards
C
  420 IFAILR = 1
      CALL X04BBF(INPUT,RREC,IFAILR)
      IF (IFAILR.NE.0) GO TO 500
C
      NLIST = NLIST + 1
C     If user has requested it, print data
      IF (MPSLST) CALL X04BAF(IPRINT,' '//RREC)
C
C     Ignore comments
      IF (RREC(1:1).EQ.'*') GO TO 420
C     check for $
      IF (RREC(15:15).EQ.DOLLAR) GO TO 420
C
      IF (RREC(1:1).NE.' ') THEN
         IF (RREC(1:6).EQ.'ENDATA') GO TO 440
C
         IF ( .NOT. MPSLST) GO TO 480
         IF (IFAIL.EQ.0) THEN
            IFAIL = 11
            WRITE (REC,FMT=99998) NLIST
         END IF
         GO TO 420
      END IF
C
      IF (RREC(4:4).NE.' ' .OR. RREC(13:14).NE.'  ' .OR. RREC(23:24)
     *    .NE.'  ' .OR. RREC(37:39).NE.'   ') THEN
         IF ( .NOT. MPSLST) GO TO 480
         IF (IFAIL.EQ.0) THEN
            IFAIL = 11
            WRITE (REC,FMT=99998) NLIST
         END IF
         WRITE (TREC,FMT=99998) NLIST
         CALL X04BAY(IPRINT,2,TREC)
         GO TO 420
      END IF
C
      IF ( .NOT. GOTNM) THEN
         NMBND = RREC(5:12)
         GOTNM = .TRUE.
      END IF
C
      IF (RREC(5:12).NE.NMBND) GO TO 420
C
      FNDNAM = .TRUE.
C     find which column name
C
      CALL H02BUW(N,NAMES,RREC(15:22),1,N,JMARK,J)
      IF (J.LT.1) THEN
         IF (IFAIL.EQ.0) THEN
            IFAIL = 9
            WRITE (REC,FMT=99983) RREC(15:22)
         END IF
         IF ( .NOT. MPSLST) RETURN
         WRITE (TREC,FMT=99983) RREC(15:22)
         CALL X04BAY(IPRINT,2,TREC)
         GO TO 420
      END IF
C
C     Convert to number
      READ (RREC,FMT=99994,ERR=480) AD(1)
C
C     select form of constraint and do the appropriate thing
C     to the bounds vectors
C
      IF (RREC(1:4).EQ.' UP ') THEN
         BU(J) = AD(1)
         IF (AD(1).EQ.0.0D0) BL(J) = -BIGBND
      ELSE IF (RREC(1:4).EQ.' UI ') THEN
         BU(J) = AD(1)
         IF (INTVAR(N).LT.1) NUMINT = NUMINT + 1
         INTVAR(J) = 1
         IF (AD(1).EQ.0.0D0) BL(J) = -BIGBND
      ELSE IF (RREC(1:4).EQ.' BV ') THEN
         BU(J) = 1.0D0
         IF (INTVAR(J).LT.1) NUMINT = NUMINT + 1
         INTVAR(J) = 1
      ELSE IF (RREC(1:4).EQ.' LO ') THEN
         BL(J) = AD(1)
      ELSE IF (RREC(1:4).EQ.' LI ') THEN
         BL(J) = AD(1)
         IF (INTVAR(J).LT.1) NUMINT = NUMINT + 1
         INTVAR(J) = 1
      ELSE IF (RREC(1:4).EQ.' FX ') THEN
         BU(J) = AD(1)
         BL(J) = AD(1)
      ELSE IF (RREC(1:4).EQ.' FR ') THEN
         BU(J) = BIGBND
         BL(J) = -BIGBND
      ELSE IF (RREC(1:4).EQ.' MI ') THEN
         BU(J) = 0.0D0
         BL(J) = -BIGBND
      ELSE IF (RREC(1:4).EQ.' PL ') THEN
         BU(J) = BIGBND
         BL(J) = 0.0D0
      ELSE
C        bound type not allowed
         IF (IFAIL.EQ.0) THEN
            IFAIL = 8
            WRITE (REC,FMT=99982) RREC(1:4)
         END IF
         IF ( .NOT. MPSLST) RETURN
         WRITE (TREC,FMT=99982) RREC(1:4)
         CALL X04BAY(IPRINT,2,TREC)
      END IF
      GO TO 420
C
C     **************************************************************
C
C     ENDATA
C     set variables
  440 DO 460 J = 1, N
         X(J) = 0.0D0
         IF (J.EQ.2*(J/2)) X(J) = 1.0D0
  460 CONTINUE
C
C     check for BOUNDS name
C
      IF ( .NOT. FNDNAM) THEN
         IF (IFAIL.EQ.0) THEN
            IFAIL = 14
            WRITE (REC,FMT=99991) NMBND
         END IF
         IF ( .NOT. MPSLST) RETURN
         WRITE (TREC,FMT=99991) NMBND
         CALL X04BAY(IPRINT,2,TREC)
      END IF
C
C     Print number of integer variables, if requested.
      IF (NUMINT.GT.0 .AND. MPSLST) THEN
         WRITE (RREC,FMT=99989) NUMINT
         CALL X04BAF(IPRINT,RREC)
      END IF
C     If user has requested it, print data
      IF (MPSLST) THEN
         WRITE (TREC,FMT=99981)
         CALL X04BAY(IPRINT,3,TREC)
         WRITE (TREC,FMT=99980) NMOBJ, NMRHS
         CALL X04BAY(IPRINT,2,TREC)
         WRITE (TREC,FMT=99979) NMRNG, NMBND
         CALL X04BAY(IPRINT,2,TREC)
      END IF
      RETURN
C
C     **************************************************
C     read error before end of file
C
  480 IF (IFAIL.EQ.0) THEN
         IFAIL = 11
         WRITE (REC,FMT=99998) NLIST
      END IF
      RETURN
C
C     **************************************************
C     end of file found before the ENDATA card
C
  500 IF (IFAIL.EQ.0) THEN
         IFAIL = 10
         WRITE (REC,FMT=99999)
      END IF
      RETURN
C
C     **************************************************
C
99999 FORMAT (' ** The last line must be the ENDATA indicator line',/)
99998 FORMAT (' ** line',I5,3X,'is not a comment nor a valid line',/)
99997 FORMAT (' ** column name with leading blank or non-alphanumeric ',
     *       'character',2X,A8,/)
99996 FORMAT (' ** integer marker is not in a correct position',/)
99995 FORMAT (' ** row name with leading blank or non-alphanumeric cha',
     *       'racter',2X,A8,/)
99994 FORMAT (24X,D12.0,13X,D12.0)
99993 FORMAT (' ** RHS name was not found',2X,A8,/)
99992 FORMAT (' ** RANGES name was not found',2X,A8,/)
99991 FORMAT (' ** BOUNDS name was not found',2X,A8,/)
99990 FORMAT (/' MPSX INPUT LISTING',/' ------------------')
99989 FORMAT (11X,'NUMBER OF INTEGER VARIABLES = ',I3)
99988 FORMAT (' ** illegal constraint type -- ( ',A4,')',/)
99987 FORMAT (' ** duplicate row name  ',A8,'  - ignored',/)
99986 FORMAT (' ** no rows specified',/)
99985 FORMAT (' ** row name  ',A8,' not defined in ROWS section',/)
99984 FORMAT (' ** no columns specified',/)
99983 FORMAT (' ** column name (',A8,') is not defined in the COLUMNS ',
     *       'section',/)
99982 FORMAT (' ** illegal bound type (',A4,')',/)
99981 FORMAT (/' NAMES SELECTED',/' --------------')
99980 FORMAT (' OBJECTIVE',6X,A8,/' RHS',12X,A8)
99979 FORMAT (' RANGES',9X,A8,/' BOUNDS',9X,A8)
99978 FORMAT (' ** no objective function found',/)
99977 FORMAT (' ** too many rows. Limit is',I5,', but the actual numbe',
     *       'r required is',I5,/)
99976 FORMAT (' ** too many columns. Limit is',I5,', but the actual nu',
     *       'mber required is',I5,/)
      END
