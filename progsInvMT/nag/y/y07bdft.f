      SUBROUTINE Y07BDF(JOB,M,N,A,LDA,TITLE,NSIG,NCOLS,NOUT)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C  Purpose
C  =======
C
C  Y07BDF prints out a matrix,  or its transpose,  using an  E format on
C  device  NOUT.  Each element is printed to  NSIG  significant figures.
C  A heading may be printed before the matrix.
C
C  Rectangular,  trapezoidal and  symmetric matrices  may be  printed by
C  this routine.
C
C  Whether the matrix,  or its transpose is printed  and whether  or not
C  elements above or below the  diagonal  are  printed  is determined by
C  the parameter JOB.
C
C  Parameters
C  ==========
C
C  JOB    - CHARACTER*1.
C
C           On entry,  JOB  specifies the  type of printing  as follows.
C
C           JOB = 'N' or 'n' or 'R' or 'r'
C                ( No transpose, Row by row )
C
C              The M by N part of A is printed row by row.
C
C           JOB = 'C' or 'c' or 'T' or 't'
C                ( Column by column, Transpose )
C
C              The M by N part of A' is printed row by row, so that A is
C              printed column by column.
C
C           JOB = 'L' or 'l'
C                ( Lower trapezoidal )
C
C              The  M by min( M, N )  lower trapezoidal  part  of  A  is
C              printed  row by row.  In this case  the  N by min( M, N )
C              strictly upper trapezoidal part of  A  is not referenced.
C
C           JOB = 'U' or 'u'
C                ( Upper trapezoidal )
C
C              The  min( M, N ) by N  upper trapezoidal  part  of  A  is
C              printed  column by column,  so that the  N by min( M, N )
C              lower trapezoidal part of  A'  is printed  row by row. In
C              this case the  M by min( M, N ) strictly lower triangular
C              part of A is not referenced.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C
C           On entry, M specifies the number of rows of A.  M must be at
C           least 1.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N specifies the number of columns of A.  N must be
C           at least 1.
C
C           Unchanged on exit.
C
C  A      - REAL             array of DIMENSION ( LDA, n ).
C
C           Before entry,  A  must contain  the  matrix  to be  printed.
C           If  M is greater than 1 and there is not room on the line to
C           print a complete row, the heading  'Row i' is printed before
C           the i(th) row of A (or A').
C
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C
C           On entry, LDA specifies the first dimension of A as declared
C           in the  calling  (sub) program.  LDA  must  be at  least  M.
C
C           Unchanged on exit.
C
C  TITLE  - CHARACTER  of  length  at least  1,  but  only  one  record.
C
C           On entry,  TITLE  specifies a  heading  to be printed before
C           printing the matrix  A.  Unless  TITLE  is blank,  TITLE  is
C           printed using the format
C
C              FORMAT( 1X, A )
C
C           and a  blank line  is  printed  following  the heading.  Any
C           trailing blanks in TITLE are not printed.
C
C           Unchanged on exit.
C
C  NSIG   - INTEGER.
C
C           On entry,  NSIG  specifies the number of significant figures
C           to which the elements of  A are printed. The field width for
C           each element of  A will be  ( NSIG + 7 ). If NSIG is outside
C           the range  ( 1, 92 )  then the  value 7  is used in place of
C           NSIG.
C
C           Unchanged on exit.
C
C  NCOLS  - INTEGER.
C
C           On entry,  NCOLS  specifies  the  maximum number of printing
C           positions per line. If NCOLS is not large enough to allow at
C           least one element per line,  or if NCOLS is greater than 132
C           then the value 72 is used in place of NCOLS.
C
C           Unchanged on exit.
C
C  NOUT   - INTEGER.
C
C           On entry,  NOUT specifies the device number for printing. If
C           NOUT  is negative then the value returned by the Nag Library
C           routine X04ABF is used in place of NOUT.
C
C           Unchanged on exit.
C
C  Further comments
C  ================
C
C  To print the elements of the 6 by 5 part of a matrix A DIMENSIONed as
C  a  10 by 8 array,  to 4 significant figures,  on unit 6, with no more
C  than  60 print positions  per line and the heading  'Matrix A',  then
C  Y07BDF may be called as
C
C     CALL Y07BDF( 'No transpose', 6, 5, A, 10, 'Matrix A', 4, 60, 6 )
C
C  To get  the default values  for  NSIG, NCOLS and NOUT,  Y07BDF may be
C  called as
C
C     CALL Y07BDF( 'No transpose', 6, 5, A, 10, 'Matrix A', -1, -1, -1 )
C
C  If A is part of a matrix B partitioned as
C
C     A = ( B1  B2 ) ,
C         ( B3  A  )
C
C  where  B1 is an  l by k  matrix ( l.ge.0, k.ge.0 ), then this routine
C  may be called with the parameter  A as  B( L + 1, K + 1 ) and  LDA as
C  the  first dimension of  B  as declared in the calling (sub) program.
C
C  Nag Fortran 77 auxilliary linear algebra routine.
C
C  -- Written on 2-December-1982.
C     Sven Hammarling, NAG Central Office.
C     This version dated 16-december-1987.
C
C
C     .. Scalar Arguments ..
      INTEGER           LDA, M, N, NCOLS, NOUT, NSIG
      CHARACTER*1       JOB
      CHARACTER*(*)     TITLE
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*)
C     .. Local Scalars ..
      INTEGER           I, II, IOUT, J, JJ, K, LENT, NC, NS
      LOGICAL           L
      CHARACTER*12      FRMAT
      CHARACTER*133     REC
C     .. External Subroutines ..
      EXTERNAL          X04BAF, Y07BAZ
C     .. Intrinsic Functions ..
      INTRINSIC         LEN, MIN
C     .. Executable Statements ..
      IOUT = NOUT
      NS = NSIG
      NC = NCOLS
C
C     Y07BAZ checks NOUT, NSIG and NCOLS and returns the format.
C     K gives the number of elements to be printed per line.
C
      CALL Y07BAZ(NS,NC,IOUT,FRMAT)
      K = NC/(NS+7)
C
C     L is true if a row or column does not fit onto a single line.
C
      L = (((JOB.EQ.'N') .OR. (JOB.EQ.'n') .OR. (JOB.EQ.'R')
     *    .OR. (JOB.EQ.'r')) .AND. (N.GT.K) .AND. (M.GT.1))
     *    .OR. (((JOB.EQ.'L') .OR. (JOB.EQ.'l')) .AND. (MIN(M,N).GT.K)
     *     .AND. (M.GT.1)) .OR. (((JOB.EQ.'C') .OR. (JOB.EQ.'c')
     *    .OR. (JOB.EQ.'T') .OR. (JOB.EQ.'t')) .AND. (M.GT.K)
     *    .AND. (N.GT.1)) .OR. (((JOB.EQ.'U') .OR. (JOB.EQ.'u'))
     *    .AND. (MIN(M,N).GT.K) .AND. (N.GT.1))
C
C     Print the heading.
C
      DO 20 LENT = LEN(TITLE), 1, -1
         IF (TITLE(LENT:LENT).NE.' ') GO TO 40
   20 CONTINUE
   40 CONTINUE
      IF (LENT.GT.0) THEN
         WRITE (REC,FMT=99999) TITLE(1:LENT)
         CALL X04BAF(IOUT,REC)
         IF ( .NOT. L) CALL X04BAF(IOUT,' ')
      END IF
C
C     Print the matrix.
C
      IF ((JOB.EQ.'N') .OR. (JOB.EQ.'n') .OR. (JOB.EQ.'R')
     *    .OR. (JOB.EQ.'r')) THEN
         DO 80 I = 1, M
            IF (L) THEN
               WRITE (REC,FMT=99998) I
               CALL X04BAF(IOUT,' ')
               CALL X04BAF(IOUT,REC)
            END IF
            DO 60 JJ = 1, N, K
               WRITE (REC,FMT=FRMAT) (A(I,J),J=JJ,MIN(JJ+K-1,N))
               CALL X04BAF(IOUT,REC)
   60       CONTINUE
   80    CONTINUE
      ELSE IF ((JOB.EQ.'L') .OR. (JOB.EQ.'l')) THEN
         DO 120 I = 1, M
            IF (L) THEN
               WRITE (REC,FMT=99998) I
               CALL X04BAF(IOUT,' ')
               CALL X04BAF(IOUT,REC)
            END IF
            DO 100 JJ = 1, MIN(I,N), K
               WRITE (REC,FMT=FRMAT) (A(I,J),J=JJ,MIN(JJ+K-1,MIN(I,N)))
               CALL X04BAF(IOUT,REC)
  100       CONTINUE
  120    CONTINUE
      ELSE IF ((JOB.EQ.'C') .OR. (JOB.EQ.'c') .OR. (JOB.EQ.'T')
     *         .OR. (JOB.EQ.'t')) THEN
         DO 160 J = 1, N
            IF (L) THEN
               WRITE (REC,FMT=99998) J
               CALL X04BAF(IOUT,' ')
               CALL X04BAF(IOUT,REC)
            END IF
            DO 140 II = 1, M, K
               WRITE (REC,FMT=FRMAT) (A(I,J),I=II,MIN(II+K-1,M))
               CALL X04BAF(IOUT,REC)
  140       CONTINUE
  160    CONTINUE
      ELSE IF ((JOB.EQ.'U') .OR. (JOB.EQ.'u')) THEN
         DO 200 J = 1, N
            IF (L) THEN
               WRITE (REC,FMT=99998) J
               CALL X04BAF(IOUT,' ')
               CALL X04BAF(IOUT,REC)
            END IF
            DO 180 II = 1, MIN(J,M), K
               WRITE (REC,FMT=FRMAT) (A(I,J),I=II,MIN(II+K-1,MIN(J,M)))
               CALL X04BAF(IOUT,REC)
  180       CONTINUE
  200    CONTINUE
      END IF
      CALL X04BAF(IOUT,' ')
      CALL X04BAF(IOUT,' ')
      RETURN
C
C     End of Y07BDF. ( SMPRE  )
C
99999 FORMAT (1X,A)
99998 FORMAT (' Row ',I5)
      END
