      SUBROUTINE Y07CAZ(NSIG,NCOLS,NOUT,FRMAT)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C  Purpose
C  =======
C
C  Y07CAZ  returns,  in the  character  FRMAT,  a  format  of  the  form
C
C     ( nn( 1P Enn.nn, 1H,, 1P Enn.nn ) )
C
C  where each n is a digit.
C
C  Parameters
C  ==========
C
C  NSIG  - INTEGER.
C
C          Before  entry,  NSIG  specifies  the  number  of  significant
C          figures required. The two digits following each decimal point
C          will be  ( NSIG - 1 )  and the  previous  two digits  will be
C          ( NSIG + 7 )  for the  first value and  ( NSIG + 6 )  for the
C          second value.
C
C          If  NSIG  is not in the range  ( 1, 92 )  NSIG  is set to  7,
C          otherwise NSIG is unchanged on exit.
C
C  NCOLS - INTEGER.
C
C          Before  entry,  NCOLS  must contain  the  number  of printing
C          positions per line.  The first two digits  in the format will
C          be such that no more than  NCOLS print positions per line are
C          used. For example, if NSIG = 4 and NCOLS = 72 then FRMAT will
C          return the format
C
C             ( 3( 1P E11.3, 1H,, 1P E10.3 ) )
C
C          If  NCOLS is not positive or if  NCOLS is such that the first
C          two digits  are  larger than  99  or are not positive,  or if
C          NCOLS  is greater than 132 then NCOLS is set to 72, otherwise
C          NCOLS is unchanged on exit.
C
C          Note that the first print position in the format is blank and
C          so the format is safe to use when the first print position is
C          used for carriage control.
C
C  NOUT  - INTEGER.
C
C          Before entry, NOUT must contain a device number for printing.
C          Note  that  this  routine  does  not  perform  any  printing.
C          If NOUT is negative then NOUT is set to the value returned by
C          the Nag Library routine  X04ABF, otherwise  NOUT is unchanged
C          on exit.
C
C  FRMAT - CHARACTER*27.
C          On return FRMAT will contain the required format in character
C          positions 1 to 27.
C
C  Nag Fortran 77 basic linear algebra routine.
C
C  -- Written on 28-April-1983.
C     Sven Hammarling.
C
C
C     .. Scalar Arguments ..
      INTEGER           NCOLS, NOUT, NSIG
      CHARACTER*27      FRMAT
C     .. Local Scalars ..
      INTEGER           I, J, M, N
C     .. External Functions ..
      CHARACTER*1       Y07BAY
      EXTERNAL          Y07BAY
C     .. External Subroutines ..
      EXTERNAL          X04ABF
C     .. Executable Statements ..
      IF (NOUT.LT.0) CALL X04ABF(0,NOUT)
      IF ((NSIG.LT.1) .OR. (NSIG.GT.92)) NSIG = 7
      IF ((NCOLS.LT.1) .OR. (NCOLS.GT.132)) NCOLS = 72
C
      FRMAT = '(  (1PE  .  ,1H,,1PE  .  ))'
      I = NSIG + 7
      M = I/10
      N = I - 10*M
      IF (M.EQ.0) M = 10
C
      FRMAT(8:8) = Y07BAY(M)
      FRMAT(9:9) = Y07BAY(N)
C
      J = I - 1
      M = J/10
      N = J - 10*M
      IF (M.EQ.0) M = 10
C
      FRMAT(21:21) = Y07BAY(M)
      FRMAT(22:22) = Y07BAY(N)
C
      J = NSIG - 1
      M = J/10
      N = J - 10*M
      IF (M.EQ.0) M = 10
C
      FRMAT(11:11) = Y07BAY(M)
      FRMAT(24:24) = FRMAT(11:11)
      FRMAT(12:12) = Y07BAY(N)
      FRMAT(25:25) = FRMAT(12:12)
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
C     End of Y07CAZ. ( KSETCE )
C
      END
