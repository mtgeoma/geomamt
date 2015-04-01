      SUBROUTINE Y07DAZ(NFIG,NCOLS,NOUT,FRMAT)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C  Purpose
C  =======
C
C  Y07DAZ  returns,  in  the  character  FRMAT,  a  format  of the  form
C
C     ( nnInn )
C
C  where each n is a digit.
C
C  Parameters
C  ==========
C
C  NFIG  - INTEGER.
C
C          Before entry,  NFIG  specifies the number of figures required
C          in  each integer.  The  two digits  following  the  I  in the
C          format  will be  ( NFIG + 3 ).
C
C          If  NFIG is not in the range ( 1, 30 ) then NFIG is set to 5.
C          Otherwise NFIG is unchanged on exit.
C
C  NCOLS - INTEGER.
C
C          Before entry,  NCOLS  must  contain  the  number  of printing
C          positions per line.  The  first two digits in the format will
C          be such that no more than  NCOLS print positions per line are
C          used. For example, if NFIG = 4 and NCOLS = 72 then FRMAT will
C          return the format
C
C             ( 10I7 )
C
C          If  NCOLS is not positive, or if NCOLS is such that the first
C          two digits  are larger than  99  or  are not positive,  or if
C          NCOLS is greater than 132 then  NCOLS is set to 72, otherwise
C          NCOLS is unchanged on exit.
C
C          Note  that the  first print position  in the  format  will be
C          blank,  so  long  as  no integer  printed  with  this  format
C          actually has more than  NFIG digits, so that in this case the
C          format  is safe to use when  the first print position is used
C          for carriage control.
C
C  NOUT  - INTEGER.
C
C          Before entry, NOUT must contain a device number for printing.
C          Note that this routine does not perform any printing. If NOUT
C          is not positive then NOUT is set to the value returned by the
C          Nag Library routine  X04ABF,  otherwise  NOUT is unchanged on
C          exit.
C
C  FRMAT - CHARACTER*7.
C          On return FRMAT will contain the required format in character
C          positions 1 to 7.
C
C  Nag Fortran 77 basic linear algebra routine.
C
C  -- Written on 8-October-1984.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Scalar Arguments ..
      INTEGER           NCOLS, NFIG, NOUT
      CHARACTER*7       FRMAT
C     .. Local Scalars ..
      INTEGER           I, J, M, N
C     .. External Functions ..
      CHARACTER*1       Y07BAY
      EXTERNAL          Y07BAY
C     .. External Subroutines ..
      EXTERNAL          X04ABF
C     .. Executable Statements ..
      IF (NOUT.LT.0) CALL X04ABF(0,NOUT)
      IF ((NFIG.LT.1) .OR. (NFIG.GT.30)) NFIG = 5
      IF ((NCOLS.LT.1) .OR. (NCOLS.GT.132)) NCOLS = 72
C
      FRMAT = '(  I  )'
      I = NFIG + 3
      M = I/10
      N = I - 10*M
      IF (M.EQ.0) M = 10
C
      FRMAT(5:5) = Y07BAY(M)
      FRMAT(6:6) = Y07BAY(N)
C
      J = NCOLS/I
      IF ((J.LT.1) .OR. (J.GT.99)) THEN
         NCOLS = 72
         J = NCOLS/I
         IF (J.LT.1) J = 1
      END IF
      M = J/10
      N = J - 10*M
      IF (M.EQ.0) M = 10
C
      FRMAT(2:2) = Y07BAY(M)
      FRMAT(3:3) = Y07BAY(N)
C
      RETURN
C
C     End of Y07DAZ. ( KSETI  )
C
      END
