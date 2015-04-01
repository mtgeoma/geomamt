      SUBROUTINE Y07BBZ(NFIGB,NFIGA,NCOLS,NOUT,FRMAT)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C  Purpose
C  =======
C
C  Y07BBZ  returns,  in the  character  FRMAT,  a  format  of  the  form
C
C     ( nnFnn.nn )
C
C  where each n is a digit.
C
C  Parameters
C  ==========
C
C  NFIGB - INTEGER.
C
C          Before entry,  NFIGB specifies the number of figures required
C          before the decimal point.
C
C          If NFIGB is not in the range ( 0, 30 ), or if NFIGB and NFIGA
C          are both zero then  NFIGB  is set to  5.  Otherwise  NFIGB is
C          unchanged on exit.
C
C  NFIGA - INTEGER.
C
C          Before entry,  NFIGA specifies the number of figures required
C          after  the decimal point.  The final two digits in the format
C          will  be   NFIGA   and   the  previous  two  digits  will  be
C          ( NFIGB + NFIGA + 3 ).
C          If  NFIGA is not in the range ( 0, 30 ) then  NFIGA is set to
C          5. Otherwise NFIGA is unchanged on exit.
C
C  NCOLS - INTEGER.
C
C          Before entry,  NCOLS  must  contain  the  number  of printing
C          positions per line.  The first two digits  in the format will
C          be such that no more than  NCOLS print positions per line are
C          used.  For  example,  if  NFIGB = 2, NFIGA = 5 and NCOLS = 72
C          then FRMAT will return the format
C
C             ( 7F10.5 )
C
C          If  NCOLS is not positive, or if NCOLS is such that the first
C          two digits  are larger than  99  or are  not positive,  or if
C          NCOLS is greater than 132 then  NCOLS is set to 72, otherwise
C          NCOLS is unchanged on exit.
C
C          Note  that  the first print position  in the  format  will be
C          blank  so long  as no  value printed out  using  this  format
C          actually has more than  ( NFIGB + NFIGA ) digits,  so that in
C          this case  the format  is safe  to use  when the  first print
C          position is used for carriage control.
C
C  NOUT  - INTEGER.
C
C          Before entry, NOUT must contain a device number for printing.
C          Note  that  this  routine  does  not  perform  any  printing.
C
C          If  NOUT  is  not positive  then  NOUT  is set  to the  value
C          returned by the Nag Library routine  X04ABF,  otherwise  NOUT
C          is unchanged on exit.
C
C  FRMAT - CHARACTER*10.
C          On return FRMAT will contain the required format in character
C          positions 1 to 10.
C
C  Nag Fortran 77 basic linear algebra routine.
C
C  -- Written on 22-February-1984.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Scalar Arguments ..
      INTEGER           NCOLS, NFIGA, NFIGB, NOUT
      CHARACTER*10      FRMAT
C     .. Local Scalars ..
      INTEGER           I, J, M, N
C     .. External Functions ..
      CHARACTER*1       Y07BAY
      EXTERNAL          Y07BAY
C     .. External Subroutines ..
      EXTERNAL          X04ABF
C     .. Executable Statements ..
      IF (NOUT.LT.0) CALL X04ABF(0,NOUT)
      IF ((NFIGB.LT.0) .OR. (NFIGB.GT.30)) NFIGB = 5
      IF ((NFIGA.LT.0) .OR. (NFIGA.GT.30)) NFIGA = 5
      IF ((NFIGB.EQ.0) .AND. (NFIGA.EQ.0)) NFIGB = 5
      IF ((NCOLS.LT.1) .OR. (NCOLS.GT.132)) NCOLS = 72
C
      FRMAT = '(  F  .  )'
      M = NFIGA/10
      N = NFIGA - 10*M
      IF (M.EQ.0) M = 10
C
      FRMAT(8:8) = Y07BAY(M)
      FRMAT(9:9) = Y07BAY(N)
C
      I = NFIGB + NFIGA + 3
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
C     End of Y07BBZ. ( KSETSF )
C
      END
