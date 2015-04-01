      SUBROUTINE G05CAY(REINIT)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     called by G05CAF, G05FAF, G05FBF or G05FDF when needed, to fill
C     the internal array RV in COMMON block CG05CA with new
C     pseudo-random numbers.
C
C     G05CAY uses a multiplicative congruential algorithm
C
C     N := N * 13**13 modulo 2**59
C
C     where N is a notional variable internal to G05CAY. The value of N
C     is converted to a real number in the range 0.0 to 1.0 by scaling
C     by 2**(-59), with care taken that the result lies strictly
C     between 0.0 and 1.0.
C
C     N is initially set to 123456789*(2**32+1) but can be changed
C     by a call to G05CBF or G05CCF.
C
C     G05CAY generates number 63 at a time, in order to achieve
C     efficiency on vector-processing machines. The first call of
C     G05CAY generates 63 consecutive values of N, N(i), i = 1,...,63.
C     Subsequent calls generate the next set of 63 values of N by
C
C     N(i) := N(i) * (13**13)**63 modulo 2**59, for i = 1,...,63.
C
C     The value 63 is defined as the symbol LV in a parameter statement
C     in each routine which needs it. The particular value 63 was
C     chosen because of special properties of the multiplier
C     (13**13)**63 modulo 2**59, which permit efficient multi-length
C     arithmetic when ILIM = 4 (see below). Only a few values of LV
C     have such properties.
C
C     ******************** ADVICE FOR IMPLEMENTORS *********************
C
C     These notes are intended to guide implementors through the text
C     changes necessary to implement the basic random number generator
C     routines G05CAY, G05CAZ, G05CBF, G05CCF, G05CFZ, G05CGZ. Please
C     follow these guidelines, and consult NAG Central Office if in any
C     doubt or difficulty. Please send a listing of your final text for
C     these routines to Central Office.
C
C     1.  Read "DETAILS-NOTE-1" below.
C         Decide the relevant value of ILIM, say nn, taking account of
C         the suggestion for 'long' integers.
C
C     2.  Activate all lines beginning CAnn.
C
C     3.  Read "DETAILS-NOTE-2" below.
C         Check whether your compiler has the functions ISHFT and IAND
C         (or equivalent functions) and compiles inline code for them.
C
C     4.  If ISHFT and IAND or equivalent functions are available as
C         inline functions, activate all lines beginning CYnn. If
C         necessary, change the function names. Otherwise activate all
C         lines beginning CXnn.
C
C     ******************************************************************
C
C     ************************ DETAILS-NOTE-1 **************************
C
C     The algorithm requires that the values of N and of the multi-
C     plier 13**13 be stored as 59-bit unsigned integers and that
C     the least significant 59 bits of their product be computed. On
C     most machines this can be done much more efficiently in
C     machine code than in Fortran. The Fortran code given here is
C     intended to give guidance on a machine code implementation,
C     and to provide a less efficient implementation as a fall-back.
C
C     The 59-bit integer N is stored as a multiple-length integer in
C     the array B. In fact for convenience the 60-bit integer 2*N is
C     stored. The multiplier 13**13 is stored in the array M.
C     The multiplier (13**13)**63 modulo 2**59 is stored in the array
C     MLV in exactly the same way as the basic multiplier is stored in
C     the array M.
C
C     The number of elements in N and M (ILIM) and the number of bits
C     used in each element of N and M (IBITS) depend on the number
C     of bits (including sign) in an integer variable as follows -
C
C        ILIM     IBITS     number of bits in integer variable
C          4        15                 .ge. 32
C          3        20                 .ge. 41
C          2        30                 .ge. 60
C
C     For greatest efficiency ILIM should be chosen as small as
C     possible.
C
C     N.B. the most significant bits of N are stored in B(I,ILIM),
C     the next most significant bits in B(I,ILIM-1), . . . , and
C     the least significant bits in B(I,1). The multiplier is stored
C     in M(ILIM), M(ILIM-1), . . . , M(1) in the same way.
C
C     Note -
C
C     1) in the above table the value of IBITS is less than half the
C     number of bits in an integer variable. This ensures that the
C     necessary integer products can be formed and summed correctly
C     without integer overflow. However many machines have instruc-
C     tions for forming double-length integer products. A machine
C     code implementation can take advantage of this and allow IBITS
C     to be as large (or almost as large) as the number of bits in
C     an integer variable and ILIM to be correspondingly smaller.
C     This should be much more efficient.
C
C     2) the figures in the rightmost column in the above table are
C     correct for the specific value of the multiplier. They are
C     certainly not correct for arbitrary 60-bit arithmetic.
C
C     3) it may well be advantageous to use 'long' integers, if
C     available, within G05CAY, even if they are not used
C     elsewhere in the library.
C
C     Variant code for the array declarations and data statements
C     is supplied in comments beginning CAnn where the digits nn are
C     the value of ILIM.
C
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LV
      PARAMETER         (LV=63)
      INTEGER           ILIM
      PARAMETER         (ILIM=4)
CA03  PARAMETER         (ILIM=3)
CA02  PARAMETER         (ILIM=2)
      DOUBLE PRECISION  ONE, R2
      PARAMETER         (ONE=1.0D0,R2=0.5D0)
      DOUBLE PRECISION  RP1, RP2
      PARAMETER         (RP1=R2**60,RP2=R2**30)
CA03  DOUBLE PRECISION  RP1, RP2
CA03  PARAMETER         (RP1=R2**60,RP2=R2**40)
CA02  DOUBLE PRECISION  RP1
CA02  PARAMETER         (RP1=R2**60)
C     .. Scalar Arguments ..
      LOGICAL           REINIT
C     .. Scalars in Common ..
      INTEGER           DEFOPT, OPTION, POSSOP, KV
C     .. Arrays in Common ..
      DOUBLE PRECISION  RV(LV)
      INTEGER           B(0:LV,ILIM)
C     .. Local Scalars ..
      DOUBLE PRECISION  ONEM
      INTEGER           I, T1, T2, T3, U, V
CX03  INTEGER           I, T1, T2, U, V
CX02  INTEGER           I, T1, U, V
CY04  INTEGER           I, T1, T2, T3, T4
CY03  INTEGER           I, T1, T2, T3
CY02  INTEGER           I, T1, T2
      LOGICAL           INIT
C     .. Local Arrays ..
      INTEGER           M(ILIM), MLV(ILIM)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. Intrinsic Functions ..
      INTRINSIC         SIGN
C     .. Common blocks ..
      COMMON            /AG05CA/B, OPTION, POSSOP, DEFOPT
      COMMON            /CG05CA/RV, KV
C     .. Save statement ..
      SAVE              /AG05CA/, /CG05CA/, ONEM, INIT
C     .. Data statements ..
      DATA              INIT / .TRUE. /
      DATA              M /
     *                  17917, 13895, 19930,     8 /
CA03 *                  247293, 485810,    275 /
CA02 *                  455329277,    282074 /
      DATA              MLV /
     *                  85,  3703,  6070,  6822 /
CA03 *                  753749, 972915, 218309 /
CA02 *                  121339989, 223549366 /
C     .. Executable Statements ..
C
C     ************************ DETAILS-NOTE-2 **************************
C
C     It is advantageous to use non-standard Fortran intrinsic
C     functions for shifting and masking if these are available and if
C     they are compiled as in-line code without the overhead of a
C     subroutine call. Alternative code is given which uses the integer
C     functions:
C
C     ISHFT(I,J) to shift I J bits to the left (a negative value of
C                 J indicating a right shift)
C     IAND(I,J)  to form the logical and of I and J
C
C     It may be necesssary to replace these by calls to different
C     intrinsic functions provided by the fortran compiler.
C
C     Variant code for this computation is supplied in comments
C     beginning CXnn (using only arithmetic operations) or in
C     comments beginning CYnn (using shifting and masking functions)
C     where the digits nn are the value of ILIM.
C
C     ******************************************************************
C
      IF (INIT.OR.REINIT) THEN
         INIT = .FALSE.
         ONEM = ONE - X02AJF()
C
C        Generate first buffer of LV integers by multiplying
C        recursively by M modulo 2**59.
C        This loop cannot be vectorized.
C
         DO 20 I = 1, LV
            V = B(I-1,1)*M(1)
            U = V/32768
            B(I,1) = V - 32768*U
            V = U + B(I-1,2)*M(1) + B(I-1,1)*M(2)
            U = V/32768
            B(I,2) = V - 32768*U
            V = U + B(I-1,3)*M(1) + B(I-1,2)*M(2) + B(I-1,1)*M(3)
            U = V/32768
            B(I,3) = V - 32768*U
            V = U + B(I-1,4)*M(1) + B(I-1,3)*M(2) + B(I-1,2)*M(3)
     *            + B(I-1,1)*M(4)
            U = V/32768
            B(I,4) = V - 32768*U
CX03        V = B(I-1,1)*M(1)
CX03        U = V/1048576
CX03        B(I,1) = V - 1048576*U
CX03        V = U + B(I-1,2)*M(1) + B(I-1,1)*M(2)
CX03        U = V/1048576
CX03        B(I,2) = V - 1048576*U
CX03        V = U + B(I-1,3)*M(1) + B(I-1,2)*M(2) + B(I-1,1)*M(3)
CX03        U = V/1048576
CX03        B(I,3) = V - 1048576*U
CX02        V = B(I-1,1)*M(1)
CX02        U = V/1073741824
CX02        B(I,1) = V - 1073741824*U
CX02        V = U + B(I-1,2)*M(1) + B(I-1,1)*M(2)
CX02        U = V/1073741824
CX02        B(I,2) = V - 1073741824*U
CY04        T1 = B(I-1,1)*M(1)
CY04        T2 = ISHFT(T1,-15) + B(I-1,2)*M(1) + B(I-1,1)*M(2)
CY04        T3 = ISHFT(T2,-15) + B(I-1,3)*M(1) + B(I-1,2)*M(2)
CY04 *                         + B(I-1,1)*M(3)
CY04        T4 = ISHFT(T3,-15) + B(I-1,4)*M(1) + B(I-1,3)*M(2)
CY04 *                         + B(I-1,2)*M(3) + B(I-1,1)*M(4)
CY04        B(I,4) = IAND(T4,32767)
CY04        B(I,3) = IAND(T3,32767)
CY04        B(I,2) = IAND(T2,32767)
CY04        B(I,1) = IAND(T1,32767)
CY03        T1 = B(I-1,1)*M(1)
CY03        T2 = ISHFT(T1,-20) + B(I-1,2)*M(1) + B(I-1,1)*M(2)
CY03        T3 = ISHFT(T2,-20) + B(I-1,3)*M(1) + B(I-1,2)*M(2)
CY03 *                         + B(I-1,1)*M(3)
CY03        B(I,3) = IAND(T3,1048575)
CY03        B(I,2) = IAND(T2,1048575)
CY03        B(I,1) = IAND(T1,1048575)
CY02        T1 = B(I-1,1)*M(1)
CY02        T2 = ISHFT(T1,-30) + B(I-1,2)*M(1) + B(I-1,1)*M(2)
CY02        B(I,2) = IAND(T2,1073741823)
CY02        B(I,1) = IAND(T1,1073741823)
   20    CONTINUE
      ELSE
C
C        Generate next buffer of LV integers by multiplying in
C        parallel by M**LV modulo 2**59.
C
         DO 40 I = 1, LV
            V = B(I,1)*MLV(1)
            U = V/32768
            T1 = V - 32768*U
            V = U + B(I,2)*MLV(1) + B(I,1)*MLV(2)
            U = V/32768
            T2 = V - 32768*U
            V = U + B(I,3)*MLV(1) + B(I,2)*MLV(2) + B(I,1)*MLV(3)
            U = V/32768
            T3 = V - 32768*U
            V = U + B(I,4)*MLV(1) + B(I,3)*MLV(2) + B(I,2)*MLV(3)
     *            + B(I,1)*MLV(4)
            U = V/32768
            B(I,4) = V - 32768*U
            B(I,3) = T3
            B(I,2) = T2
            B(I,1) = T1
CX03        V = B(I,1)*MLV(1)
CX03        U = V/1048576
CX03        T1 = V - 1048576*U
CX03        V = U + B(I,2)*MLV(1) + B(I,1)*MLV(2)
CX03        U = V/1048576
CX03        T2 = V - 1048576*U
CX03        V = U + B(I,3)*MLV(1) + B(I,2)*MLV(2) + B(I,1)*MLV(3)
CX03        U = V/1048576
CX03        B(I,3) = V - 1048576*U
CX03        B(I,2) = T2
CX03        B(I,1) = T1
CX02        V = B(I,1)*MLV(1)
CX02        U = V/1073741824
CX02        T1 = V - 1073741824*U
CX02        V = U + B(I,2)*MLV(1) + B(I,1)*MLV(2)
CX02        U = V/1073741824
CX02        B(I,2) = V - 1073741824*U
CX02        B(I,1) = T1
CY04        T1 = B(I,1)*MLV(1)
CY04        T2 = ISHFT(T1,-15) + B(I,2)*MLV(1) + B(I,1)*MLV(2)
CY04        T3 = ISHFT(T2,-15) + B(I,3)*MLV(1) + B(I,2)*MLV(2)
CY04 *                         + B(I,1)*MLV(3)
CY04        T4 = ISHFT(T3,-15) + B(I,4)*MLV(1) + B(I,3)*MLV(2)
CY04 *                         + B(I,2)*MLV(3) + B(I,1)*MLV(4)
CY04        B(I,4) = IAND(T4,32767)
CY04        B(I,3) = IAND(T3,32767)
CY04        B(I,2) = IAND(T2,32767)
CY04        B(I,1) = IAND(T1,32767)
CY03        T1 = B(I,1)*MLV(1)
CY03        T2 = ISHFT(T1,-20) + B(I,2)*MLV(1) + B(I,1)*MLV(2)
CY03        T3 = ISHFT(T2,-20) + B(I,3)*MLV(1) + B(I,2)*MLV(2)
CY03 *                         + B(I,1)*MLV(3)
CY03        B(I,3) = IAND(T3,1048575)
CY03        B(I,2) = IAND(T2,1048575)
CY03        B(I,1) = IAND(T1,1048575)
CY02        T1 = B(I,1)*MLV(1)
CY02        T2 = ISHFT(T1,-30) + B(I,2)*MLV(1) + B(I,1)*MLV(2)
CY02        B(I,2) = IAND(T2,1073741823)
CY02        B(I,1) = IAND(T1,1073741823)
   40    CONTINUE
      END IF
C
C     Convert integers in B to real numbers in (0.0,1.0) stored in RV.
C
      DO 60 I = 1, LV
         RV(I) = MIN(ONEM,(B(I,4)*32768+B(I,3))*RP2
     *                     +(B(I,2)*32768+B(I,1))*RP1)
CX03     RV(I) = MIN(ONEM,(B(I,3)*1048576+B(I,2))*RP2 + B(I,1)*RP1)
CX02     RV(I) = MIN(ONEM,(B(I,2)*1073741824+B(I,1))*RP1)
CY04     RV(I) = MIN(ONEM,(ISHFT(B(I,4),15)+B(I,3))*RP2
CY04 *                     +(ISHFT(B(I,2),15)+B(I,1))*RP1)
CY03     RV(I) = MIN(ONEM,(ISHFT(B(I,3),20)+B(I,2))*RP2 + B(I,1)*RP1)
CY02     RV(I) = MIN(ONEM,(ISHFT(B(I,2),30)+B(I,1))*RP1)
   60 CONTINUE
      KV = 0
C
      RETURN
      END
