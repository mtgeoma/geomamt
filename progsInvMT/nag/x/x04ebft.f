      SUBROUTINE X04EBF(MATRIX,DIAG,M,N,A,LDA,FORMAT,TITLE,LABROW,RLABS,
     *                  LABCOL,CLABS,NCOLS,INDENT,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 16A REVISED. IER-1050 (JUN 1993).
C     Prints a general integer matrix.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='X04EBF')
      INTEGER           ONE, ZERO
      PARAMETER         (ONE=1,ZERO=0)
C     .. Scalar Arguments ..
      INTEGER           IFAIL, INDENT, LDA, M, N, NCOLS
      CHARACTER*1       DIAG, LABCOL, LABROW, MATRIX
      CHARACTER*(*)     FORMAT, TITLE
C     .. Array Arguments ..
      INTEGER           A(LDA,*)
      CHARACTER*(*)     CLABS(*), RLABS(*)
C     .. Local Scalars ..
      INTEGER           AA, CL, CLBWID, CLEFT, CRIGHT, CSHIFT, DIAGEL,
     *                  FINIS2, FINISH, I, IERR, INCOLS, INDNT, J, K,
     *                  LWID, MAXEL, MINEL, ND, NELEMS, NELS, NOUT,
     *                  NREC, NSETS, NTITLE, NUMWID, OFFSET, RLBWID,
     *                  START, START2
      LOGICAL           GENERL, LOWER, PRDIAG
      CHARACTER*10      ICFORM, IRFORM
      CHARACTER*81      FORMT
      CHARACTER*87      FORM
      CHARACTER*132     BLANKS, TILTER
      CHARACTER*133     INFILE
      CHARACTER*266     INFIL2
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          X04ABF, X04BAF, X04CBZ
C     .. Intrinsic Functions ..
      INTRINSIC         LEN, MAX, MIN
C     .. Executable Statements ..
C
      IERR = 0
      GENERL = MATRIX .EQ. 'G' .OR. MATRIX .EQ. 'g'
      LOWER = MATRIX .EQ. 'L' .OR. MATRIX .EQ. 'l'
      IF (NCOLS.LE.0 .OR. NCOLS.GT.132) THEN
         INCOLS = 80
      ELSE
         INCOLS = NCOLS
      END IF
      IF (INDENT.LT.0 .OR. INDENT.GE.INCOLS) THEN
         INDNT = 0
      ELSE
         INDNT = INDENT
      END IF
      INCOLS = INCOLS - INDNT
      BLANKS = ' '
C
C     Check for incorrect arguments.
      IF ( .NOT. GENERL .AND. .NOT. LOWER .AND. MATRIX.NE.'U' .AND.
     *    MATRIX.NE.'u') THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99999) MATRIX
      ELSE IF ((LABROW.NE.'N' .AND. LABROW.NE.'n' .AND. LABROW.NE.
     *         'I' .AND. LABROW.NE.'i' .AND. LABROW.NE.'C' .AND.
     *         LABROW.NE.'c') .OR. (LABCOL.NE.'N' .AND. LABCOL.NE.
     *         'n' .AND. LABCOL.NE.'I' .AND. LABCOL.NE.'i' .AND.
     *         LABCOL.NE.'C' .AND. LABCOL.NE.'c')) THEN
         IERR = 6
         NREC = 2
         WRITE (REC,FMT=99996) LABROW, LABCOL
      ELSE IF (M.GT.LDA) THEN
         IERR = 3
         NREC = 1
         WRITE (REC,FMT=99998) M, LDA
      ELSE IF ( .NOT. GENERL) THEN
         IF (DIAG.NE.'U' .AND. DIAG.NE.'u' .AND. DIAG.NE.'N' .AND.
     *       DIAG.NE.'n' .AND. DIAG.NE.'B' .AND. DIAG.NE.'b') THEN
            IERR = 2
            NREC = 1
            WRITE (REC,FMT=99997) MATRIX, DIAG
         END IF
      END IF
      IF (IERR.NE.0) GO TO 260
C
C     Get the advisory message unit number.
      CALL X04ABF(0,NOUT)
C
      FORMT = FORMAT
      IF (FORMT.EQ.' ') THEN
C        Construct the FORMAT code of smallest field width that is
C        large enough to print all matrix elements.
C        First find the largest and smallest elements to be printed.
         MAXEL = ZERO
         MINEL = ZERO
         IF (GENERL) THEN
            DO 40 J = 1, N
               DO 20 I = 1, M
                  AA = A(I,J)
                  IF (AA.GT.MAXEL) MAXEL = AA
                  IF (AA.LT.MINEL) MINEL = AA
   20          CONTINUE
   40       CONTINUE
         ELSE
            IF (DIAG.EQ.'N' .OR. DIAG.EQ.'n') THEN
               DO 50 J = 1, MIN(M,N)
                  AA = A(J,J)
                  IF (AA.GT.MAXEL) MAXEL = AA
                  IF (AA.LT.MINEL) MINEL = AA
   50          CONTINUE
            END IF
            DO 100 J = 1, N
               IF (LOWER) THEN
                  DO 60 I = J + 1, M
                     AA = A(I,J)
                     IF (AA.GT.MAXEL) MAXEL = AA
                     IF (AA.LT.MINEL) MINEL = AA
   60             CONTINUE
               ELSE
                  DO 80 I = 1, J - 1
                     AA = A(I,J)
                     IF (AA.GT.MAXEL) MAXEL = AA
                     IF (AA.LT.MINEL) MINEL = AA
   80             CONTINUE
               END IF
  100       CONTINUE
         END IF
C
         WRITE (INFILE,FMT='(I40)') MAXEL
         CALL X04CBZ(INFILE,START,FINISH)
         WRITE (INFILE,FMT='(I40)') MINEL
         CALL X04CBZ(INFILE,START2,FINIS2)
         K = MAX(FINISH-START,FINIS2-START2) + 2
         FORMT = 'I  '
         WRITE (FORMT(2:),FMT='(I2)') K
      END IF
C
C     Construct the format statement to be used internally.
      CALL X04CBZ(FORMT,START,FINISH)
      IF (FINISH-START+1.GT.80) THEN
C        The length of FORMAT is too great.
         IERR = 4
         NREC = 2
         WRITE (REC,FMT=99995) FORMAT(START:START+74)
         GO TO 260
      END IF
C
      FORM = '(999('//FORMT(START:FINISH)//'))'
C
C     Decide how wide each column of numbers is going to be,
C     by writing the number 0 to an internal file and measuring
C     the width needed. Since the width may include trailing blanks,
C     we also write 0 twice with the same format, and compute
C     the required field width from the two widths. Note that if
C     FORMAT has more than one edit descriptor in it, with different
C     field widths, for example FORMAT = 'I12,I13', then it is
C     not possible to compute a sensible value, in which case
C     the columns of the output matrix will be a little skew.
      WRITE (INFILE,FMT=FORM,ERR=120) ZERO
      CALL X04CBZ(INFILE,START,FINISH)
      WRITE (INFIL2,FMT=FORM,ERR=120) ZERO, ZERO
      CALL X04CBZ(INFIL2,START,FINIS2)
      NUMWID = FINIS2 - FINISH
C     NUMWID is the width of a number as printed using FORMT.
      GO TO 140
  120 CONTINUE
C     The format in FORM caused an error when used to print a number.
      IERR = 5
      NREC = 2
      WRITE (REC,FMT=99994) FORMAT(1:MIN(75,LEN(FORMAT)))
      GO TO 260
  140 CONTINUE
C
C     What kind of row labelling is required?
      IF (LABROW.EQ.'N' .OR. LABROW.EQ.'n') THEN
C        No row labelling.
         RLBWID = 1
      ELSE IF (LABROW.EQ.'I' .OR. LABROW.EQ.'i') THEN
C        Numeric row labelling.
         WRITE (INFILE,FMT='(I16)') M
         CALL X04CBZ(INFILE,START,FINISH)
         RLBWID = FINISH - START + 2
         IRFORM = '(I    )'
         WRITE (IRFORM(3:6),FMT='(I4)') RLBWID
      ELSE
C        User supplied row labelling.
         RLBWID = 1
         DO 160 I = 1, M
            CALL X04CBZ(RLABS(I),START,FINISH)
            RLBWID = MAX(RLBWID,FINISH-START+2)
  160    CONTINUE
      END IF
C
C     What kind of column labelling is required?
      IF (LABCOL.EQ.'I' .OR. LABCOL.EQ.'i') THEN
C        Numeric column labelling.
         WRITE (INFILE,FMT='(I16)') N
         CALL X04CBZ(INFILE,START,FINISH)
         CLBWID = FINISH - START + 2
         ICFORM = '(999I    )'
         WRITE (ICFORM(6:9),FMT='(I4)') NUMWID
      ELSE IF (LABCOL.NE.'N' .AND. LABCOL.NE.'n') THEN
C        User supplied column labelling.
         CLBWID = LEN(CLABS(1))
      END IF
C
      NELEMS = (INCOLS-1-RLBWID)/NUMWID
      IF (NELEMS.LT.1) THEN
         IERR = 7
         NREC = 2
         WRITE (REC,FMT=99993) INCOLS + INDNT, INDNT
         GO TO 260
      END IF
C     NELEMS is the number of elements that can fit into INCOLS columns.
C
      NSETS = (N-1)/NELEMS + 1
C     NSETS is the number of pieces that the matrix must be split into.
C
C     Print the title, splitting it up if more than INCOLS-1 characters.
      CALL X04CBZ(TITLE,START,FINISH)
      IF (FINISH.NE.0) THEN
         NTITLE = (FINISH-1)/(INCOLS-1) + 1
         DO 180 I = 1, NTITLE - 1
            TILTER = BLANKS(1:INDNT+1)//TITLE((I-1)*(INCOLS-1)
     *               +1:I*(INCOLS-1))
            CALL X04BAF(NOUT,TILTER)
  180    CONTINUE
         TILTER = BLANKS(1:INDNT+1)//TITLE((NTITLE-1)*(INCOLS-1)
     *            +1:FINISH)
         CALL X04BAF(NOUT,TILTER)
      END IF
C
C     Exit after printing the title if M or N is less than 1.
      IF (M.LT.1 .OR. N.LT.1) GO TO 260
C
C     Print the matrix, with row and column labels if requested.
      CSHIFT = 0
C     CSHIFT is the offset into the current set of columns, when
C     the matrix cannot be printed in one go but has to be split.
      DO 240 I = 1, NSETS
         IF (I.EQ.NSETS) THEN
            NELS = N - (NSETS-1)*NELEMS
         ELSE
            NELS = NELEMS
         END IF
         IF (LABCOL.EQ.'I' .OR. LABCOL.EQ.'i') THEN
C           Construct the numeric column labels.
            INFILE = ' '
            WRITE (INFILE(RLBWID+2+INDNT:),FMT=ICFORM) (CL,CL=CSHIFT+1,
     *        CSHIFT+NELS)
         ELSE IF (LABCOL.NE.'N' .AND. LABCOL.NE.'n') THEN
C           Process the user-supplied column labels.
            INFILE = ' '
            LWID = MIN(CLBWID,NUMWID)
            OFFSET = RLBWID + 1 + INDNT
            DO 200 K = CSHIFT + 1, CSHIFT + NELS
               CALL X04CBZ(CLABS(K)(1:LWID),START,FINISH)
               IF (START.EQ.0) THEN
                  START = 1
                  FINISH = 1
               END IF
               INFILE(OFFSET+NUMWID-FINISH+START:
     *           OFFSET+NUMWID+FINISH-START) = CLABS(K) (START:FINISH)
               OFFSET = OFFSET + NUMWID
  200       CONTINUE
         END IF
C        Output the column labels.
         IF (LABCOL.NE.'N' .AND. LABCOL.NE.'n') THEN
            CALL X04BAF(NOUT,INFILE(1:NCOLS))
         END IF
C
C        Now print each row in turn.
         DO 220 J = 1, M
            INFILE = ' '
C
C           Insert the row label.
            IF (LABROW.EQ.'I' .OR. LABROW.EQ.'i') THEN
               WRITE (INFILE(INDNT+1:INDNT+RLBWID),FMT=IRFORM) J
            ELSE IF (LABROW.NE.'N' .AND. LABROW.NE.'n') THEN
               CALL X04CBZ(RLABS(J),START,FINISH)
               IF (START.EQ.0) THEN
                  START = 1
                  FINISH = 1
               END IF
               INFILE(INDNT+RLBWID-FINISH+START:INDNT+RLBWID) = RLABS(J)
     *           (START:FINISH)
            END IF
C
            IF (GENERL) THEN
C              General rectangular matrix.
               WRITE (INFILE(INDNT+RLBWID+2:),FMT=FORM) (A(J,CL),
     *           CL=CSHIFT+1,CSHIFT+NELS)
            ELSE
C              Upper or lower triangular matrix.
               ND = MAX(0,J-(I-1)*NELEMS)
C              ND is the position of the Jth row diagonal element.
               IF (LOWER) THEN
                  CLEFT = CSHIFT + 1
                  CRIGHT = CSHIFT + MIN(ND-1,NELS)
               ELSE
                  CLEFT = CSHIFT + ND + 1
                  CRIGHT = CSHIFT + NELS
               END IF
C              CLEFT and CRIGHT are the leftmost and rightmost elements
C              of the current row to be printed, excluding the diagonal.
               PRDIAG = DIAG .NE. 'B' .AND. DIAG .NE. 'b' .AND. ND .GT.
     *                  0 .AND. ND .LE. NELS
C              PRDIAG is true if a diagonal element appears in the
C              current matrix row section, and it is to be printed.
               IF (PRDIAG) THEN
                  IF (DIAG.EQ.'U' .OR. DIAG.EQ.'u') THEN
                     DIAGEL = ONE
                  ELSE
                     DIAGEL = A(J,J)
                  END IF
               END IF
C
               IF (LOWER) THEN
C                 Lower triangular matrix.
                  IF (PRDIAG) THEN
                     WRITE (INFILE(INDNT+RLBWID+2:),FMT=FORM) (A(J,CL),
     *                 CL=CLEFT,CRIGHT), DIAGEL
                  ELSE
                     WRITE (INFILE(INDNT+RLBWID+2:),FMT=FORM) (A(J,CL),
     *                 CL=CLEFT,CRIGHT)
                  END IF
               ELSE
C                 Upper triangular matrix.
                  IF (PRDIAG) THEN
                     WRITE (INFILE(INDNT+RLBWID+2+NUMWID*(ND-1):),
     *                 FMT=FORM) DIAGEL, (A(J,CL),CL=CLEFT,CRIGHT)
                  ELSE
                     IF (CLEFT.LE.CRIGHT) THEN
C                       Have to do the check on CLEFT and CRIGHT to
C                       avoid INDNT+RLBWID+2+NUMWID*ND possibly being
C                       out of range.
                        WRITE (INFILE(INDNT+RLBWID+2+NUMWID*ND:),
     *                    FMT=FORM) (A(J,CL),CL=CLEFT,CRIGHT)
                     END IF
                  END IF
               END IF
            END IF
C
C           Output the (partial) matrix row.
            CALL X04BAF(NOUT,INFILE(1:NCOLS))
C
  220    CONTINUE
C
         CSHIFT = CSHIFT + NELEMS
         IF (I.NE.NSETS) THEN
            CALL X04BAF(NOUT,' ')
         END IF
  240 CONTINUE
C
  260 IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (1X,'** On entry, MATRIX is not valid : MATRIX = ''',A,
     *       '''         ')
99998 FORMAT (1X,'** On entry, M.lt.LDA : M = ',I16,', LDA = ',I16)
99997 FORMAT (1X,'** On entry, MATRIX = ''',A,''', but DIAG is not val',
     *       'id : DIAG = ''',A,'''')
99996 FORMAT (1X,'** On entry, either LABROW or LABCOL is not valid',
     *       /4X,'LABROW = ''',A,''', LABCOL = ''',A,'''.')
99995 FORMAT (1X,'** On entry, FORMAT has more than 80 characters: the',
     *       ' first 75 are',/4X,A)
99994 FORMAT (1X,'** The format specifier in FORMAT cannot be used to ',
     *       'print a number: FORMAT =',/4X,A)
99993 FORMAT (1X,'** On entry, NCOLS-INDENT is not wide enough to hold',
     *       ' at least one matrix',/4X,'column: NCOLS = ',I16,', INDE',
     *       'NT =',I16)
      END
