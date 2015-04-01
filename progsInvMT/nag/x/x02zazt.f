      SUBROUTINE X02ZAZ
C     MARK 16 REVISED. IER-1046 (JUN 1993).
C
C***********************************************************************
C
C     NAG version of the Stanford routine MCHPAR.
C     Sven Hammarling, NAG Central Office.
C
C     X02ZAZ sets machine parameters as follows:
C
C     WMACH(  1 ) = nbase  = base of floating-point arithmetic.
C     WMACH(  2 ) = ndigit = no. of base ( nbase ) digits in the
C                            mantissa
C     WMACH(  3 ) = eps    = relative machine accuracy. (X02AJF.)
C     WMACH(  4 ) = rteps  = sqrt( eps ).
C     WMACH(  5 ) = rmin   = small positive floating-point number whose
C                             reciprocal does not overflow.
C     WMACH(  6 ) = rtrmin = sqrt( rmin ).
C     WMACH(  7 ) = rmax   = 1/rmin
C     WMACH(  8 ) = rtrmax = sqrt( rmax ).
C     WMACH(  9 ) = undflw = 0 if underflow is not fatal, +ve otherwise.
C     WMACH( 10 ) = nin    = input  stream unit number. ( 5.)
C     WMACH( 11 ) = nout   = output stream unit number.
C                          = advisory message unit number. ( X04ABF.)
C     WMACH( 12 ) = nerr   = error    message unit number. ( X04AAF.)
C     WMACH( 13 )
C     WMACH( 14 )   Not currently used.
C     WMACH( 15 )
C
C     Note that constants that represent integers may hold a number just
C     less than the integer, so that the integer should be recovered by
C     adding, say, 0.25. e.g.
C
C     IBASE = WMACH( 1 ) + 0.25
C
C***********************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER        (ZERO=0.0D+0)
C     .. Arrays in Common ..
      DOUBLE PRECISION WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION EPS, RMAX, RMIN, UNDFLW
      INTEGER          NBASE, NDIGIT, NERR, NOUT
      LOGICAL          FIRST
C     .. External Functions ..
      DOUBLE PRECISION X02AJF, X02AMF
      INTEGER          X02BHF, X02BJF
      LOGICAL          X02DAF
      EXTERNAL         X02AJF, X02AMF, X02BHF, X02BJF, X02DAF
C     .. External Subroutines ..
      EXTERNAL         X04AAF, X04ABF
C     .. Intrinsic Functions ..
      INTRINSIC        SQRT
C     .. Common blocks ..
      COMMON           /AX02ZA/WMACH
C     .. Save statement ..
      SAVE             /AX02ZA/, FIRST
C     .. Data statements ..
      DATA             FIRST/.TRUE./
C     .. Executable Statements ..
C
      IF (FIRST) THEN
         FIRST = .FALSE.
C
         IF (X02DAF(ZERO)) THEN
            UNDFLW = 1
         ELSE
            UNDFLW = 0
         END IF
         NBASE = X02BHF()
         NDIGIT = X02BJF()
         EPS = X02AJF()
         RMIN = X02AMF()
         RMAX = 1/RMIN
C
         WMACH(1) = NBASE
         WMACH(2) = NDIGIT
         WMACH(3) = EPS
         WMACH(4) = SQRT(EPS)
         WMACH(5) = RMIN
         WMACH(6) = SQRT(RMIN)
         WMACH(7) = RMAX
         WMACH(8) = SQRT(RMAX)
         WMACH(9) = UNDFLW
      END IF
      CALL X04ABF(0,NOUT)
      WMACH(10) = 5
      WMACH(11) = NOUT
      CALL X04AAF(0,NERR)
      WMACH(12) = NERR
      RETURN
C
C     End of  X02ZAZ. (MCHPAR)
C
      END
