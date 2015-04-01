      SUBROUTINE Y07BAZ(NSIG,NCOLS,NOUT,FRMAT)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C  Purpose
C  =======
C
C  Y07BAZ  returns,  in  the  character  FRMAT,  a  format  of the  form
C
C     ( 1P nnEnn.nn )
C
C  where each n is a digit.
C
C  Parameters
C  ==========
C
C  NSIG  - INTEGER.
C
C          Before  entry,  NSIG  specifies  the  number  of  significant
C          figures required. The final two digits in the format will  be
C          ( NSIG - 1 )   and   the   previous   two   digits   will  be
C          ( NSIG + 7 ).
C
C          If  NSIG  is not in the range  ( 1, 92 ),  NSIG  is set to 7,
C          otherwise NSIG is unchanged on exit.
C
C  NCOLS - INTEGER.
C
C          Before entry,  NCOLS  must contain  the  number  of  printing
C          positions per line.  The first  two digits in the format will
C          be such that no more than  NCOLS print positions per line are
C          used. For example, if NSIG = 8 and NCOLS = 72 then FRMAT will
C          return the format
C
C             ( 1P 4E15.7 )
C
C          If  NCOLS is not positive, or if NCOLS is such that the first
C          two digits  are  larger than  99  or  are not positive, or if
C          NCOLS is greater than 132 then  NCOLS is set to 72, otherwise
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
C
C          If NOUT is negative then NOUT is set to the value returned by
C          the Nag Library routine  X04ABF, otherwise  NOUT is unchanged
C          on exit.
C
C  FRMAT - CHARACTER*12.
C          On return FRMAT will contain the required format in character
C          positions 1 to 12.
C
C  Nag Fortran 77 basic linear algebra routine.
C
C  -- Written on 1-December-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Scalar Arguments ..
      INTEGER           NCOLS, NOUT, NSIG
      CHARACTER*12      FRMAT
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
      FRMAT = '(1P  E  .  )'
      I = NSIG + 7
      M = I/10
      N = I - 10*M
      IF (M.EQ.0) M = 10
C
      FRMAT(7:7) = Y07BAY(M)
      FRMAT(8:8) = Y07BAY(N)
C
      J = NSIG - 1
      M = J/10
      N = J - 10*M
      IF (M.EQ.0) M = 10
C
      FRMAT(10:10) = Y07BAY(M)
      FRMAT(11:11) = Y07BAY(N)
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
      FRMAT(4:4) = Y07BAY(M)
      FRMAT(5:5) = Y07BAY(N)
C
      RETURN
C
C     End of Y07BAZ. ( KSETSE )
C
      END
