      SUBROUTINE X04DBF(MATRIX,DIAG,M,N,A,LDA,USEFRM,FORMAT,TITLE,
     *                  LABROW,RLABS,LABCOL,CLABS,NCOLS,INDENT,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 16A REVISED. IER-1049 (JUN 1993).
C     Prints a general complex matrix.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='X04DBF')
      COMPLEX*16        CZERO
      PARAMETER         (CZERO=(0.0D+0,0.0D+0))
      DOUBLE PRECISION  FRNGLW, FRNGUP, ONE, ZERO
      PARAMETER         (FRNGLW=1.0D-3,FRNGUP=9.99999991D+3,ONE=1.0D+0,
     *                  ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           IFAIL, INDENT, LDA, M, N, NCOLS
      CHARACTER         DIAG, LABCOL, LABROW, MATRIX, USEFRM
      CHARACTER*(*)     FORMAT, TITLE
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*)
      CHARACTER*(*)     CLABS(*), RLABS(*)
C     .. Local Scalars ..
      COMPLEX*16        DIAGEL
      DOUBLE PRECISION  AA, BASE, LBASE, MAXEL, MINEL, RNDIGS
      INTEGER           CL, CLBWID, CLEFT, CRIGHT, CSHIFT, FINIS2,
     *                  FINISH, I, IERR, INCOLS, INDNT, J, K, LWID, ND,
     *                  NDIGS1, NDIGS2, NELEMS, NELS, NOUT, NREC, NSETS,
     *                  NTITLE, NUMWID, OFFSET, RLBWID, START
      LOGICAL           ABOVE, BRACK, GENERL, LOWER, PRDIAG
      CHARACTER*10      ICFORM, IRFORM
      CHARACTER*81      FORMT
      CHARACTER*132     BLANKS, TILTER
      CHARACTER*133     INFBOT, INFILE
      CHARACTER*189     FORM
      CHARACTER*266     INFIL2
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF, X02BHF, X02BJF, X02BKF, X02BLF
      EXTERNAL          P01ABF, X02BHF, X02BJF, X02BKF, X02BLF
C     .. External Subroutines ..
      EXTERNAL          X04ABF, X04BAF, X04CBZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, LEN, LOG10, MAX, MIN, DBLE
C     .. Executable Statements ..
C
      IERR = 0
      GENERL = MATRIX .EQ. 'G' .OR. MATRIX .EQ. 'g'
      LOWER = MATRIX .EQ. 'L' .OR. MATRIX .EQ. 'l'
      ABOVE = USEFRM .EQ. 'A' .OR. USEFRM .EQ. 'a'
      BRACK = USEFRM .EQ. 'B' .OR. USEFRM .EQ. 'b'
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
         IERR = 7
         NREC = 2
         WRITE (REC,FMT=99996) LABROW, LABCOL
      ELSE IF ( .NOT. ABOVE .AND. .NOT. BRACK .AND. USEFRM.NE.'D' .AND.
     *         USEFRM.NE.'d') THEN
         IERR = 4
         NREC = 1
         WRITE (REC,FMT=99995) USEFRM
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
      IF (IERR.NE.0) GO TO 560
C
C     Get the advisory message unit number.
      CALL X04ABF(0,NOUT)
C
      FORMT = FORMAT
      IF (FORMT.EQ.'*') THEN
C        Construct an E FORMAT that is wide enough to distinguish
C        between adjacent machine numbers.
         BASE = X02BHF()
         LBASE = LOG10(BASE)
         RNDIGS = X02BJF()*LBASE
         NDIGS1 = RNDIGS
         IF (NDIGS1.LT.RNDIGS) NDIGS1 = NDIGS1 + 1
C        NDIGS1 is the number of significant decimal digits required
C        for the mantissa.
         RNDIGS = LOG10(MAX(-X02BKF(),X02BLF())*LBASE) + 2
         NDIGS2 = RNDIGS
         IF (NDIGS2.LT.RNDIGS) NDIGS2 = NDIGS2 + 1
C        NDIGS2 is the number of decimal places required for the
C        exponent, including letter 'E' and sign.
         FORMT = '1P,E   .   '
         WRITE (FORMT(5:7),FMT='(I3)') NDIGS1 + NDIGS2 + 3
         WRITE (FORMT(9:11),FMT='(I3)') NDIGS1 - 1
      ELSE IF (FORMT.EQ.' ') THEN
C        Construct either a fixed point FORMAT, if the elements to be
C        printed are a reasonable size, i.e. all lie inside the range
C        0.001 to 9999.9999, or a floating-point FORMAT otherwise,
C        printing to 5 significant digits.
C        First find the largest and smallest elements to be printed,
C        ignoring zeros.
         MAXEL = ONE
         MINEL = ONE
         IF (GENERL) THEN
            DO 40 J = 1, N
               DO 20 I = 1, M
                  AA = MAX(ABS(DBLE(A(I,J))),ABS(DIMAG(A(I,J))))
                  IF (AA.GT.MAXEL) MAXEL = AA
                  IF (AA.NE.ZERO .AND. AA.LT.MINEL) MINEL = AA
   20          CONTINUE
   40       CONTINUE
         ELSE
            IF (DIAG.EQ.'N' .OR. DIAG.EQ.'n') THEN
               DO 50 J = 1, MIN(M,N)
                  AA = MAX(ABS(DBLE(A(J,J))),ABS(DIMAG(A(J,J))))
                  IF (AA.GT.MAXEL) MAXEL = AA
                  IF (AA.NE.ZERO .AND. AA.LT.MINEL) MINEL = AA
   50          CONTINUE
            END IF
            DO 100 J = 1, N
               IF (LOWER) THEN
                  DO 60 I = J + 1, M
                     AA = MAX(ABS(DBLE(A(I,J))),ABS(DIMAG(A(I,J))))
                     IF (AA.GT.MAXEL) MAXEL = AA
                     IF (AA.NE.ZERO .AND. AA.LT.MINEL) MINEL = AA
   60             CONTINUE
               ELSE
                  DO 80 I = 1, J - 1
                     AA = MAX(ABS(DBLE(A(I,J))),ABS(DIMAG(A(I,J))))
                     IF (AA.GT.MAXEL) MAXEL = AA
                     IF (AA.NE.ZERO .AND. AA.LT.MINEL) MINEL = AA
   80             CONTINUE
               END IF
  100       CONTINUE
         END IF
C
         IF (MINEL.GE.FRNGLW .AND. MAXEL.LT.FRNGUP) THEN
C           If all elements to be printed are moderately sized,
C           use a fixed point FORMAT ...
            IF (MAXEL.LE.ONE) THEN
               FORMT = 'F8.4'
            ELSE
               FORMT = 'F11.4'
            END IF
         ELSE
C           ... otherwise use a floating-point FORMAT.
            FORMT = '1PE13.4'
         END IF
      END IF
C
C     Construct the format statement to be used internally.
      CALL X04CBZ(FORMT,START,FINISH)
      IF (FINISH-START+1.GT.80) THEN
C        The length of FORMAT is too great.
         IERR = 5
         NREC = 2
         WRITE (REC,FMT=99994) FORMAT(START:START+74)
         GO TO 560
      END IF
C
      IF (BRACK) THEN
C        Insert brackets and comma into the edit descriptor.
         FORM = '(999(:,'' ('','//FORMT(START:FINISH)
     *          //','','','//FORMT(START:FINISH)//','')''))'
      ELSE
         FORM = '(999(:,'//FORMT(START:FINISH)//'))'
      END IF
C
C     Decide how wide each column of numbers is going to be,
C     by writing the number ZERO to an internal file and measuring
C     the width needed. Since the width may include trailing blanks,
C     we also write ZERO twice with the same format, and compute
C     the required field width from the two widths. Note that if
C     FORMAT has more than one edit descriptors in it, with different
C     field widths, for example FORMAT = 'E12.3,E13.4', then it may
C     not be possible to compute a sensible value, in which case
C     the columns of the output matrix will be a little skew.
      IF (ABOVE) THEN
         WRITE (INFILE,FMT=FORM,ERR=120) ZERO
         WRITE (INFIL2,FMT=FORM,ERR=120) ZERO, ZERO
      ELSE
         WRITE (INFILE,FMT=FORM,ERR=120) CZERO
         WRITE (INFIL2,FMT=FORM,ERR=120) CZERO, CZERO
      END IF
      CALL X04CBZ(INFILE,START,FINISH)
      CALL X04CBZ(INFIL2,START,FINIS2)
      NUMWID = FINIS2 - FINISH
C     NUMWID is the width of a number as printed using FORMT.
      GO TO 140
  120 CONTINUE
C     The format in FORM caused an error when used to print a number.
      IERR = 6
      NREC = 2
      WRITE (REC,FMT=99993) FORMAT(1:MIN(75,LEN(FORMAT)))
      GO TO 560
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
         IERR = 8
         NREC = 2
         WRITE (REC,FMT=99992) INCOLS + INDNT, INDNT
         GO TO 560
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
      IF (M.LT.1 .OR. N.LT.1) GO TO 560
C
C     Print the matrix, with row and column labels if requested.
      CSHIFT = 0
C     CSHIFT is the offset into the current set of columns, when
C     the matrix cannot be printed in one go but has to be split.
      DO 540 I = 1, NSETS
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
         DO 520 J = 1, M
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
               IF (ABOVE) THEN
                  WRITE (INFILE(INDNT+RLBWID+2:),FMT=FORM,ERR=220)
     *              (DBLE(A(J,CL)),CL=CSHIFT+1,CSHIFT+NELS)
  220             CONTINUE
                  INFBOT(1:INDNT+RLBWID+1) = ' '
                  WRITE (INFBOT(INDNT+RLBWID+2:),FMT=FORM,ERR=240)
     *              (DIMAG(A(J,CL)),CL=CSHIFT+1,CSHIFT+NELS)
  240             CONTINUE
               ELSE
                  WRITE (INFILE(INDNT+RLBWID+2:),FMT=FORM,ERR=260)
     *              (A(J,CL),CL=CSHIFT+1,CSHIFT+NELS)
  260             CONTINUE
               END IF
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
                     IF (ABOVE) THEN
                        WRITE (INFILE(INDNT+RLBWID+2:),FMT=FORM,ERR=280)
     *                    (DBLE(A(J,CL)),CL=CLEFT,CRIGHT), DBLE(DIAGEL)
  280                   CONTINUE
                        INFBOT(1:INDNT+RLBWID+1) = ' '
                        WRITE (INFBOT(INDNT+RLBWID+2:),FMT=FORM,ERR=300)
     *                    (DIMAG(A(J,CL)),CL=CLEFT,CRIGHT),
     *                    DIMAG(DIAGEL)
  300                   CONTINUE
                     ELSE
                        WRITE (INFILE(INDNT+RLBWID+2:),FMT=FORM,ERR=320)
     *                    (A(J,CL),CL=CLEFT,CRIGHT), DIAGEL
  320                   CONTINUE
                     END IF
                  ELSE
                     IF (ABOVE) THEN
                        WRITE (INFILE(INDNT+RLBWID+2:),FMT=FORM,ERR=340)
     *                    (DBLE(A(J,CL)),CL=CLEFT,CRIGHT)
  340                   CONTINUE
                        INFBOT(1:INDNT+RLBWID+1) = ' '
                        WRITE (INFBOT(INDNT+RLBWID+2:),FMT=FORM,ERR=360)
     *                    (DIMAG(A(J,CL)),CL=CLEFT,CRIGHT)
  360                   CONTINUE
                     ELSE
                        WRITE (INFILE(INDNT+RLBWID+2:),FMT=FORM,ERR=380)
     *                    (A(J,CL),CL=CLEFT,CRIGHT)
  380                   CONTINUE
                     END IF
                  END IF
               ELSE
C                 Upper triangular matrix.
                  IF (PRDIAG) THEN
                     IF (ABOVE) THEN
                        WRITE (INFILE(INDNT+RLBWID+2+NUMWID*(ND-1):),
     *                    FMT=FORM,ERR=400) DBLE(DIAGEL),
     *                    (DBLE(A(J,CL)),CL=CLEFT,CRIGHT)
  400                   CONTINUE
                        INFBOT = ' '
                        WRITE (INFBOT(INDNT+RLBWID+2+NUMWID*(ND-1):),
     *                    FMT=FORM,ERR=420) DIMAG(DIAGEL),
     *                    (DIMAG(A(J,CL)),CL=CLEFT,CRIGHT)
  420                   CONTINUE
                     ELSE
                        WRITE (INFILE(INDNT+RLBWID+2+NUMWID*(ND-1):),
     *                    FMT=FORM,ERR=440) DIAGEL,
     *                    (A(J,CL),CL=CLEFT,CRIGHT)
  440                   CONTINUE
                     END IF
                  ELSE
                     IF (CLEFT.LE.CRIGHT) THEN
C                       Have to do the check on CLEFT and CRIGHT to
C                       avoid INDNT+RLBWID+2+NUMWID*ND possibly being
C                       out of range.
                        IF (ABOVE) THEN
                           WRITE (INFILE(INDNT+RLBWID+2+NUMWID*ND:),
     *                       FMT=FORM,ERR=460) (DBLE(A(J,CL)),CL=CLEFT,
     *                       CRIGHT)
  460                      CONTINUE
                           INFBOT = ' '
                           WRITE (INFBOT(INDNT+RLBWID+2+NUMWID*ND:),
     *                       FMT=FORM,ERR=480) (DIMAG(A(J,CL)),CL=CLEFT,
     *                       CRIGHT)
  480                      CONTINUE
                        ELSE
                           WRITE (INFILE(INDNT+RLBWID+2+NUMWID*ND:),
     *                       FMT=FORM,ERR=500) (A(J,CL),CL=CLEFT,CRIGHT)
  500                      CONTINUE
                        END IF
                     END IF
                  END IF
               END IF
            END IF
C
C           Output the (partial) matrix row.
            CALL X04BAF(NOUT,INFILE(1:NCOLS))
            IF (ABOVE) THEN
               CALL X04BAF(NOUT,INFBOT(1:NCOLS))
               IF (J.LT.M) CALL X04BAF(NOUT,' ')
            END IF
C
  520    CONTINUE
C
         CSHIFT = CSHIFT + NELEMS
         IF (I.NE.NSETS) THEN
            CALL X04BAF(NOUT,' ')
         END IF
  540 CONTINUE
C
  560 IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (1X,'** On entry, MATRIX is not valid : MATRIX = ''',A,
     *       '''')
99998 FORMAT (1X,'** On entry, M.gt.LDA : M = ',I16,', LDA = ',I16)
99997 FORMAT (1X,'** On entry, MATRIX = ''',A,''', but DIAG is not val',
     *       'id : DIAG = ''',A,'''')
99996 FORMAT (1X,'** On entry, either LABROW or LABCOL is not valid',
     *       /4X,'LABROW = ''',A,''', LABCOL = ''',A,'''.')
99995 FORMAT (1X,'** On entry, USEFRM is not valid : USEFRM = ''',A,
     *       '''         .')
99994 FORMAT (1X,'** On entry, FORMAT has more than 80 characters: the',
     *       ' first 75 are',/4X,A)
99993 FORMAT (1X,'** The format specifier in FORMAT cannot be used to ',
     *       'print a number: FORMAT =',/4X,A)
99992 FORMAT (1X,'** On entry, NCOLS-INDENT is not wide enough to hold',
     *       ' at least one matrix',/4X,'column: NCOLS = ',I16,', INDE',
     *       'NT =',I16)
      END
