      SUBROUTINE Y07DAF(N,K,INCK,TITLE,NFIG,NCOLS,NOUT)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C  Purpose
C  =======
C
C  Y07DAF prints out the n element vector k, using an I format on device
C  NOUT. A heading may be printed before the vector.
C
C  Parameters
C  ==========
C
C  N     - INTEGER.
C
C          On entry,  N  specifies  the  number of elements of  k  to be
C          printed.  If  N  is not positive  then an immediate return is
C          made.
C
C          Unchanged on exit.
C
C  K     - INTEGER  array of  DIMENSION at least ( max( 1, lx ) ), where
C          lk = 1 + ( n - 1 )*abs( INCK ).
C
C          Before entry, the incremented array K must contain the vector
C          to be printed. If INCK is positive then the elements  K( 1 ),
C          K( 1 + INCK ), ..., K( 1 + ( N - 1 )*INCK )  are printed.  If
C          INCK  is negative then the elements  K( 1 - ( N - 1 )*INCK ),
C          ..., K( 1 - INCK ), K( 1 )  are printed.
C
C          Unchanged on exit.
C
C  INCK  - INTEGER.
C
C          On entry,  INCK specifies the increment for printing. If INCK
C          is zero then an immediate return is made.
C
C          Unchanged on exit.
C
C  TITLE - CHARACTER  of  length  at  least  1,  but  only  one  record.
C
C          On entry,  TITLE  specifies  a heading  to be  printed before
C          printing  the vector  K.  Unless  TITLE  is blank,  TITLE  is
C          printed using the format
C
C             FORMAT( 1X, A )
C
C          and  a  blank line  is  printed  following  the heading.  Any
C          trailing blanks in TITLE are not printed.
C
C          Unchanged on exit.
C
C  NFIG  - INTEGER.
C
C          On entry,  NFIG  specifies the  number of figures required in
C          each element of k. The field width for each element of k will
C          be ( NFIG + 3 ).  If NFIG is not in the range ( 1, 30 ), then
C          the value 5 is used in place of NFIG.
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
C  To print the elements of the 20 element vector  k, with a format that
C  allows  4 figures  for each element on  unit 6,  with no more than 60
C  print positions  per line and the heading  'Vector k',  Y07DAF may be
C  called as
C
C     CALL Y07DAF( 20, K, 1, 'Vector k', 4, 60, 6 )
C
C  To get the  default values  for  NFIG, NCOLS and NOUT,  Y07DAF may be
C  called as
C
C     CALL Y07DAF( 20, K, 1, 'Vector k', -1, -1, -1 )
C
C  Nag Fortran 77 basic linear algebra routine.
C
C  -- Written on 8-October-1984.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Scalar Arguments ..
      INTEGER           INCK, N, NCOLS, NFIG, NOUT
      CHARACTER*(*)     TITLE
C     .. Array Arguments ..
      INTEGER           K(*)
C     .. Local Scalars ..
      INTEGER           I, IOUT, J, L, LENT, NC, NF
      CHARACTER*7       FRMAT
      CHARACTER*133     REC
C     .. External Subroutines ..
      EXTERNAL          X04BAF, Y07DAZ
C     .. Intrinsic Functions ..
      INTRINSIC         LEN, MAX, MIN
C     .. Executable Statements ..
      IF ((N.LT.1) .OR. (INCK.EQ.0)) RETURN
C
      IOUT = NOUT
      NF = NFIG
      NC = NCOLS
C
C     Y07DAZ checks NOUT, NFIG and NCOLS and returns the format.
C     L gives the number of elements to be printed per line.
C
      CALL Y07DAZ(NF,NC,IOUT,FRMAT)
      L = NC/(NF+3)
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
      IF (INCK.GT.0) THEN
         DO 60 J = 1, 1 + (N-1)*INCK, L*INCK
            WRITE (REC,FMT=FRMAT) (K(I),I=J,MIN(J+(L-1)*INCK,1+(N-1)
     *        *INCK),INCK)
            CALL X04BAF(IOUT,REC)
   60    CONTINUE
      ELSE
         DO 80 J = 1 - (N-1)*INCK, 1, L*INCK
            WRITE (REC,FMT=FRMAT) (K(I),I=J,MAX(J+(L-1)*INCK,1),INCK)
            CALL X04BAF(IOUT,REC)
   80    CONTINUE
      END IF
C
      CALL X04BAF(IOUT,' ')
      CALL X04BAF(IOUT,' ')
      RETURN
C
C     End of Y07DAF. ( IVPRI  )
C
99999 FORMAT (1X,A)
      END
