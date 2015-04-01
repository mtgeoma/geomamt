      SUBROUTINE Y07CAF(N,X,INCX,TITLE,NSIG,NCOLS,NOUT)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C  Purpose
C  =======
C
C  Y07CAF prints out the n element vector x, using an E format on device
C  NOUT.  Each element is printed to NSIG significant figures. A heading
C  may be printed before the vector.
C
C  Parameters
C  ==========
C
C  N     - INTEGER.
C
C          On entry,  N  specifies  the  number of elements of  x  to be
C          printed.  If  N  is not positive  then an immediate return is
C          made.
C
C          Unchanged on exit.
C
C  X     - COMPLEX        array of  DIMENSION at least ( max( 1, lx ) ),
C          where  lx = 1 + ( n - 1 )*abs( INCX ).
C
C          Before entry, the incremented array X must contain the vector
C          to be printed. If INCX is positive then the elements  X( 1 ),
C          X( 1 + INCX ), ..., X( 1 + ( N - 1 )*INCX )  are printed.  If
C          INCX  is negative then the elements  X( 1 - ( N - 1 )*INCX ),
C          ..., X( 1 - INCX ), X( 1 )  are printed.
C
C          Unchanged on exit.
C
C  INCX  - INTEGER.
C
C          On entry,  INCX specifies the increment for printing. If INCX
C          is zero then an immediate return is made.
C
C          Unchanged on exit.
C
C  TITLE - CHARACTER  of  length  at  least  1,  but  only  one  record.
C
C          On entry,  TITLE  specifies  a heading  to be  printed before
C          printing  the vector  x.  Unless  TITLE  is blank,  TITLE  is
C          printed using the format
C
C             FORMAT( 1X, A )
C
C          and  a  blank line  is  printed  following  the heading.  Any
C          trailing blanks in TITLE are not printed.
C
C          Unchanged on exit.
C
C  NSIG  - INTEGER.
C
C          On entry, NSIG specifies the number of significant figures to
C          which the elements of x are printed. The field width for each
C          element of  x  will be  2*( NSIG + 7 ). A comma separates the
C          real  and  imaginary parts.  If  NSIG  is  outside  the range
C          ( 1, 92 )  then  the  value  7  is used  in  place  of  NSIG.
C
C          Unchanged on exit.
C
C  NCOLS - INTEGER.
C
C          On entry,  NCOLS  specifies  the  maximum number  of printing
C          positions per line.  If NCOLS is not large enough to allow at
C          least one element per line or if  NCOLS  is greater than  132
C          then the value 72 is used in place of NCOLS.
C
C          Unchanged on exit.
C
C  NOUT  - INTEGER.
C
C          On entry,  NOUT  specifies  the  device number  for printing.
C          If  NOUT  is  negative  then the  value returned  by the  Nag
C          Library routine X04ABF is used in place of NOUT.
C
C          Unchanged on exit.
C
C  Further comments
C  ================
C
C  To print the elements of the  20 element vector  x,  to 4 significant
C  figures on unit 6,  with no more than 60 print positions per line and
C  the heading  'Vector x',  then  Y07CAF may be called as
C
C     CALL Y07CAF( 20, X, 1, 'Vector x', 4, 60, 6 )
C
C  To get the  default values for  NSIG, NCOLS and NOUT,  Y07CAF  may be
C  called as
C
C     CALL Y07CAF( 20, X, 1, 'Vector x', -1, -1, -1 )
C
C  Nag Fortran 77 basic linear algebra routine.
C
C  -- Written on 28-April-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Scalar Arguments ..
      INTEGER           INCX, N, NCOLS, NOUT, NSIG
      CHARACTER*(*)     TITLE
C     .. Array Arguments ..
      COMPLEX*16        X(*)
C     .. Local Scalars ..
      INTEGER           I, IOUT, J, K, LENT, NC, NS
      CHARACTER*27      FRMAT
      CHARACTER*133     REC
C     .. External Subroutines ..
      EXTERNAL          X04BAF, Y07CAZ
C     .. Intrinsic Functions ..
      INTRINSIC         LEN, MAX, MIN
C     .. Executable Statements ..
      IF ((N.LT.1) .OR. (INCX.EQ.0)) RETURN
C
      IOUT = NOUT
      NS = NSIG
      NC = NCOLS
C
C     Y07CAZ checks NOUT, NSIG and NCOLS and returns the format.
C     K gives the number of elements to be printed per line.
C
      CALL Y07CAZ(NS,NC,IOUT,FRMAT)
      K = NC/(2*(NS+7))
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
         CALL X04BAF(IOUT,' ')
      END IF
C
C     Print the vector.
C
      IF (INCX.GT.0) THEN
         DO 60 J = 1, 1 + (N-1)*INCX, K*INCX
            WRITE (REC,FMT=FRMAT) (X(I),I=J,MIN(J+(K-1)*INCX,1+(N-1)
     *        *INCX),INCX)
            CALL X04BAF(IOUT,REC)
   60    CONTINUE
      ELSE
         DO 80 J = 1 - (N-1)*INCX, 1, K*INCX
            WRITE (REC,FMT=FRMAT) (X(I),I=J,MAX(J+(K-1)*INCX,1),INCX)
            CALL X04BAF(IOUT,REC)
   80    CONTINUE
      END IF
C
      CALL X04BAF(IOUT,' ')
      CALL X04BAF(IOUT,' ')
      RETURN
C
C     End of Y07CAF. ( CVPRE  )
C
99999 FORMAT (1X,A)
      END
