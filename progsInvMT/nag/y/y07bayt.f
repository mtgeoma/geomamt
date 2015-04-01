      CHARACTER*1 FUNCTION Y07BAY(I)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C  Y07BAY  returns a single character from the  Fortran 77 character set
C  as follows.
C
C     I  Y07BAY     I  Y07BAY     I  Y07BAY     I  Y07BAY     I  Y07BAY
C   -12    $
C   -11    '
C   -10    :
C    -9    .        1    1       11    A       21    K       31    U
C    -8    ,        2    2       12    B       22    L       32    V
C    -7    )        3    3       13    C       23    M       33    W
C    -6    (        4    4       14    D       24    N       34    X
C    -5    /        5    5       15    E       25    O       35    Y
C    -4    *        6    6       16    F       26    P       36    Z
C    -3    -        7    7       17    G       27    Q    OTHER
C    -2    +        8    8       18    H       28    R
C    -1    =        9    9       19    I       29    S
C     0    0       10            20    J       30    T
C
C  I is not altered by this routine.
C
C  Nag Fortran 77 basic linear algebra routine.
C
C  -- Written on 2-December-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Scalar Arguments ..
      INTEGER                     I
C     .. Local Scalars ..
      INTEGER                     J
      CHARACTER*49                K
C     .. Data statements ..
      DATA                        K(1:12)/'$'':.,)(/*-+='/
      DATA                        K(13:23)/'0123456789 '/
      DATA                        K(24:36)/'ABCDEFGHIJKLM'/
      DATA                        K(37:49)/'NOPQRSTUVWXYZ'/
C     .. Executable Statements ..
      J = I + 13
      IF ((J.LT.1) .OR. (J.GT.49)) J = 23
C
      Y07BAY = K(J:J)
      RETURN
C
C     End of Y07BAY. ( KCHAR  )
C
      END
