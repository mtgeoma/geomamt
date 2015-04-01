      SUBROUTINE Y07CBZ(NFIGB,NFIGA,NCOLS,NOUT,FRMAT)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C  Purpose
C  =======
C
C  Y07CBZ  returns,  in the  character  FRMAT,  a  format  of  the  form
C
C     ( nn( Fnn.nn, 1H,, Fnn.nn ) )
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
C          after the  decimal point.  The final two digits in the format
C          will  be   NFIGA   and  the  previous   two  digits  will  be
C          (  NFIGB  +  NFIGA  +  3  )   for   the   first   value   and
C          (  NFIGB  +  NFIGA  +  2  )   for   the   second  value.
C
C          If  NFIGA is not in the range ( 0, 30 ) then  NFIGA is set to
C          5. Otherwise NFIGA is unchanged on exit.
C
C  NCOLS - INTEGER.
C
C          Before entry,  NCOLS  must contain  the  number  of  printing
C          positions per line.  The first two digits  in the format will
C          be such that no more than  NCOLS print positions per line are
C          used.  For  example,  if  NFIGB = 2, NFIGA = 3 and NCOLS = 72
C          then FRMAT will return the format
C
C             ( 4( F8.3, 1H,, F7.3 ) )
C
C          If  NCOLS is not positive, or if NCOLS is such that the first
C          two digits  are larger than  99  or are  not positive,  or if
C          NCOLS is greater than 132 then  NCOLS is set to 72, otherwise
C          NCOLS is unchanged on exit.
C
C          Note  that  the  first print position  in the  format will be
C          blank  so long  as  no  real value  printed  by  this  format
C          actuallly has more than  ( NFIGB + NFIGA ) digits, so that in
C          this case  the  format  is  safe to use when  the first print
C          position is used for carriage control.
C
C  NOUT  - INTEGER.
C
C          Before entry, NOUT must contain a device number for printing.
C          Note that this routine does not perform any printing. If NOUT
C          is not positive then NOUT is set to the value returned by the
C          Nag Library routine  X04ABF,  otherwise  NOUT is unchanged on
C          exit.
C
C  FRMAT - CHARACTER*23.
C          On return FRMAT will contain the required format in character
C          positions 1 to 23.
C
C  Nag Fortran 77 basic linear algebra routine.
C
C  -- Written on 22-February-1984.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Scalar Arguments ..
      INTEGER           NCOLS, NFIGA, NFIGB, NOUT
      CHARACTER*23      FRMAT
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
      FRMAT = '(  (F  .  ,1H,,F  .  ))'
      M = NFIGA/10
      N = NFIGA - 10*M
      IF (M.EQ.0) M = 10
C
      FRMAT(9:9) = Y07BAY(M)
      FRMAT(20:20) = FRMAT(9:9)
      FRMAT(10:10) = Y07BAY(N)
      FRMAT(21:21) = FRMAT(10:10)
C
      I = NFIGB + NFIGA + 2
      M = I/10
      N = I - 10*M
      IF (M.EQ.0) M = 10
C
      FRMAT(17:17) = Y07BAY(M)
      FRMAT(18:18) = Y07BAY(N)
C
      I = NFIGB + NFIGA + 3
      M = I/10
      N = I - 10*M
      IF (M.EQ.0) M = 10
C
      FRMAT(6:6) = Y07BAY(M)
      FRMAT(7:7) = Y07BAY(N)
C
      J = NCOLS/(2*I)
      IF ((J.LT.1) .OR. (J.GT.99)) THEN
         NCOLS = 72
         J = NCOLS/(2*I)
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
C     End of Y07CBZ. ( KSETCF )
C
      END
