      DOUBLE PRECISION FUNCTION S17AFF(X,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 7C REVISED IER-185 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-758 (DEC 1989).
C     BESSEL FUNCTION J1(X)
C
C     **************************************************************
C
C     TO EXTRACT THE CORRECT CODE FOR A PARTICULAR MACHINE-RANGE,
C     ACTIVATE THE STATEMENTS CONTAINED IN COMMENTS BEGINNING  CDD ,
C     WHERE  DD  IS THE APPROXIMATE NUMBER OF SIGNIFICANT DECIMAL
C     DIGITS REPRESENTED BY THE MACHINE
C     DELETE THE ILLEGAL DUMMY STATEMENTS OF THE FORM
C     * EXPANSION (NNNN) *
C
C     NOTE - THE VALUE OF THE CONSTANT XBIG (DEFINED IN A DATA
C     STATEMENT) SHOULD BE REDUCED IF NECESSARY TO ENSURE THAT SIN
C     AND COS WILL RETURN A RESULT WITHOUT AN EXECUTION ERROR FOR
C     ALL ARGUMENTS X SUCH THAT ABS(X) .LE. XBIG
C
C     **************************************************************
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='S17AFF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 A, B, C, CX, G, SX, T, T2, TBPI,
     *                                 XBIG, XVSMAL, Y
C     .. Local Arrays ..
      CHARACTER*1                      P01REC(1)
C     .. External Functions ..
      INTEGER                          P01ABF
      EXTERNAL                         P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, SIGN, COS, SIN, SQRT
C     .. Data statements ..
C08   DATA XBIG,XVSMAL,TBPI/1.0D+7,1.0D-4,6.36619772D-1/
C09   DATA XBIG,XVSMAL,TBPI/1.0D+8,3.2D-5,6.366197724D-1/
C12   DATA XBIG,XVSMAL,TBPI/1.0D+11,1.0D-6,6.366197723676D-1/
C15   DATA XBIG,XVSMAL,TBPI/1.0D+14,3.2D-8,6.366197723675813D-1/
      DATA XBIG,XVSMAL,TBPI/1.0D+16,3.2D-9,6.36619772367581343D-1/
C19   DATA XBIG,XVSMAL,TBPI/1.0D+18,3.2D-10,6.3661977236758134308D-1/
C     .. Executable Statements ..
C
      T = ABS(X)
C     ERROR 1 TEST
      IF (T.GT.XBIG) GO TO 60
      IFAIL = 0
C     X RANGE TEST
      IF (T.GT.8.0D0) GO TO 40
C     SMALL X
      Y = 4.0D0
C      TEST FOR VERY SMALL X
      IF (T.LE.XVSMAL) GO TO 20
      T = 3.125D-2*T*T - 1.0D0
      T2 = 2.0D0*T
C
C      * EXPANSION (0031) *
C
C     EXPANSION (0031) EVALUATED AS Y(T)  --PRECISION 08E.09
C08   A = +2.94970701D-8
C08   B = T2*A - 7.61758781D-7
C08   C = T2*B - A + 1.58870192D-5
C08   A = T2*C - B - 2.60444389D-4
C08   B = T2*A - C + 3.24027018D-3
C08   C = T2*B - A - 2.91755248D-2
C08   A = T2*C - B + 1.77709117D-1
C08   B = T2*A - C - 6.61443934D-1
C08   C = T2*B - A + 1.28799410D+0
C08   A = T2*C - B - 1.19180116D+0
C08   Y = T*A - C + 6.48358771D-1
C
C     EXPANSION (0031) EVALUATED AS Y(T)  --PRECISION 09E.10
C09   A = -9.424212982D-10
C09   B = T2*A + 2.949707007D-8
C09   C = T2*B - A - 7.617587805D-7
C09   A = T2*C - B + 1.588701924D-5
C09   B = T2*A - C - 2.604443893D-4
C09   C = T2*B - A + 3.240270183D-3
C09   A = T2*C - B - 2.917552481D-2
C09   B = T2*A - C + 1.777091172D-1
C09   C = T2*B - A - 6.614439341D-1
C09   A = T2*C - B + 1.287994099D+0
C09   B = T2*A - C - 1.191801161D+0
C09   Y = T*B - A + 6.483587706D-1
C
C     EXPANSION (0031) EVALUATED AS Y(T)  --PRECISION 12E.13
C12   A = -5.777404200000D-13
C12   B = T2*A + 2.528123664000D-11
C12   C = T2*B - A - 9.424212981600D-10
C12   A = T2*C - B + 2.949707007278D-8
C12   B = T2*A - C - 7.617587805400D-7
C12   C = T2*B - A + 1.588701923993D-5
C12   A = T2*C - B - 2.604443893486D-4
C12   B = T2*A - C + 3.240270182684D-3
C12   C = T2*B - A - 2.917552480615D-2
C12   A = T2*C - B + 1.777091172397D-1
C12   B = T2*A - C - 6.614439341345D-1
C12   C = T2*B - A + 1.287994098858D+0
C12   A = T2*C - B - 1.191801160541D+0
C12   Y = T*A - C + 6.483587706053D-1
C
C     EXPANSION (0031) EVALUATED AS Y(T)  --PRECISION 15E.16
C15   A = -1.955400000000000D-16
C15   B = T2*A + 1.138572000000000D-14
C15   C = T2*B - A - 5.777404200000000D-13
C15   A = T2*C - B + 2.528123664000000D-11
C15   B = T2*A - C - 9.424212981600000D-10
C15   C = T2*B - A + 2.949707007278000D-8
C15   A = T2*C - B - 7.617587805400300D-7
C15   B = T2*A - C + 1.588701923993213D-5
C15   C = T2*B - A - 2.604443893485807D-4
C15   A = T2*C - B + 3.240270182683857D-3
C15   B = T2*A - C - 2.917552480615421D-2
C15   C = T2*B - A + 1.777091172397283D-1
C15   A = T2*C - B - 6.614439341345433D-1
C15   B = T2*A - C + 1.287994098857678D+0
C15   C = T2*B - A - 1.191801160541217D+0
C15   Y = T*C - B + 6.483587706052649D-1
C
C     EXPANSION (0031) EVALUATED AS Y(T)  --PRECISION 17E.18
      A = +2.95000000000000000D-18
      B = T2*A - 1.95540000000000000D-16
      C = T2*B - A + 1.13857200000000000D-14
      A = T2*C - B - 5.77740420000000000D-13
      B = T2*A - C + 2.52812366400000000D-11
      C = T2*B - A - 9.42421298160000000D-10
      A = T2*C - B + 2.94970700727800000D-8
      B = T2*A - C - 7.61758780540030000D-7
      C = T2*B - A + 1.58870192399321300D-5
      A = T2*C - B - 2.60444389348580680D-4
      B = T2*A - C + 3.24027018268385747D-3
      C = T2*B - A - 2.91755248061542077D-2
      A = T2*C - B + 1.77709117239728283D-1
      B = T2*A - C - 6.61443934134543253D-1
      C = T2*B - A + 1.28799409885767762D+0
      A = T2*C - B - 1.19180116054121687D+0
      Y = T*A - C + 6.48358770605264921D-1
C
C     EXPANSION (0031) EVALUATED AS Y(T)  --PRECISION 19E.20
C19   A = -4.0000000000000000000D-20
C19   B = T2*A + 2.9500000000000000000D-18
C19   C = T2*B - A - 1.9554000000000000000D-16
C19   A = T2*C - B + 1.1385720000000000000D-14
C19   B = T2*A - C - 5.7774042000000000000D-13
C19   C = T2*B - A + 2.5281236640000000000D-11
C19   A = T2*C - B - 9.4242129816000000000D-10
C19   B = T2*A - C + 2.9497070072780000000D-8
C19   C = T2*B - A - 7.6175878054003000000D-7
C19   A = T2*C - B + 1.5887019239932130000D-5
C19   B = T2*A - C - 2.6044438934858068000D-4
C19   C = T2*B - A + 3.2402701826838574700D-3
C19   A = T2*C - B - 2.9175524806154207660D-2
C19   B = T2*A - C + 1.7770911723972828328D-1
C19   C = T2*B - A - 6.6144393413454325277D-1
C19   A = T2*C - B + 1.2879940988576776204D+0
C19   B = T2*A - C - 1.1918011605412168725D+0
C19   Y = T*B - A + 6.4835877060526492084D-1
C
   20 S17AFF = Y*X*0.125D0
      GO TO 80
C
C     LARGE X
   40 G = T - 1.5D0/TBPI
      Y = SIGN(SQRT(TBPI/T),X)
      CX = COS(G)*Y
      SX = -SIN(G)*Y*8.0D0/T
      T = 128.0D0/(T*T) - 1.0D0
C
C      * EXPANSION (0033) *
C
C     EXPANSION (0033) EVALUATED AS Y(T)  --PRECISION 08E.09
C08   Y = ((((-1.49751260D-8)*T+2.47105358D-7)*T-7.95959347D-6)
C08  *    *T+8.98804504D-4)*T + 1.00090703D+0
C
C     EXPANSION (0033) EVALUATED AS Y(T)  --PRECISION 09E.10
C09   Y = ((((-1.497512599D-8)*T+2.471053584D-7)*T-7.959593475D-6)
C09  *    *T+8.988045041D-4)*T + 1.000907026D+0
C
C     EXPANSION (0033) EVALUATED AS Y(T)  --PRECISION 12E.13
C12   Y = (((((((+3.007485120000D-11)*T-1.825556364800D-10)
C12  *    *T+1.358072796000D-9)*T-1.470129253816D-8)
C12  *    *T+2.453682941886D-7)*T-7.959696162530D-6)
C12  *    *T+8.988049416226D-4)*T + 1.000907026278D+0
C
C     EXPANSION (0033) EVALUATED AS Y(T)  --PRECISION 15E.16
C15   Y = (((((((((((+1.092403200000000D-13)*T-3.697254400000000D-13)
C15  *    *T+1.095472640000000D-12)*T-5.071493120000000D-12)
C15  *    *T+2.723452416000000D-11)*T-1.713727974400000D-10)
C15  *    *T+1.360296919680000D-9)*T-1.470849844856000D-8)
C15  *    *T+2.453676633377600D-7)*T-7.959694699684780D-6)
C15  *    *T+8.988049416705182D-4)*T + 1.000907026278082D+0
C
C     EXPANSION (0033) EVALUATED AS Y(T)  --PRECISION 17E.18
      Y = (((((((((((((+1.24928000000000000D-14)
     *    *T-3.54508800000000000D-14)*T+6.86387200000000000D-14)
     *    *T-2.63372800000000000D-13)*T+1.14622464000000000D-12)
     *    *T-5.19113984000000000D-12)*T+2.72040729600000000D-11)
     *    *T-1.71310758400000000D-10)*T+1.36030580128000000D-9)
     *    *T-1.47085129889600000D-8)*T+2.45367662227560000D-7)
     *    *T-7.95969469843846000D-6)*T+8.98804941670557880D-4)*T +
     *    1.00090702627808217D+0
C
C     EXPANSION (0033) EVALUATED AS Y(T)  --PRECISION 19E.20
C19   Y = ((((((((((((((((-6.5536000000000000000D-16)
C19  *    *T+1.9660800000000000000D-15)*T-2.1299200000000000000D-15)
C19  *    *T+5.1200000000000000000D-15)*T-2.3080960000000000000D-14)
C19  *    *T+7.9697920000000000000D-14)*T-2.8263424000000000000D-13)
C19  *    *T+1.1377766400000000000D-12)*T-5.1772390400000000000D-12)
C19  *    *T+2.7207528960000000000D-11)*T-1.7131578496000000000D-10)
C19  *    *T+1.3603050755200000000D-9)*T-1.4708512133280000000D-8)
C19  *    *T+2.4536766229476000000D-7)*T-7.9596946984927400000D-6)
C19  *    *T+8.9880494167055608000D-4)*T + 1.0009070262780821665D+0
C
C      * EXPANSION (0034) *
C
C     EXPANSION (0034) EVALUATED AS G(T)  --PRECISION 08E.09
C08   G = (((((-7.49818190D-10)*T+6.58335466D-9)*T-8.29018528D-8)
C08  *    *T+1.82113970D-6)*T-9.62145905D-5)*T + 4.67768740D-2
C
C     EXPANSION (0034) EVALUATED AS G(T)  --PRECISION 09E.10
C09   G = (((((-7.498181901D-10)*T+6.583354662D-9)*T-8.290185280D-8)
C09  *    *T+1.821139697D-6)*T-9.621459047D-5)*T + 4.677687403D-2
C
C     EXPANSION (0034) EVALUATED AS G(T)  --PRECISION 12E.13
C12   G = (((((((((-1.167680000000D-12)*T+4.603874560000D-12)
C12  *    *T-1.826434048000D-11)*T+1.032792521600D-10)
C12  *    *T-7.152283142400D-10)*T+6.420379003440D-9)
C12  *    *T-8.291958561664D-8)*T+1.821201819899D-6)
C12  *    *T-9.621458822163D-5)*T + 4.677687402745D-2
C
C     EXPANSION (0034) EVALUATED AS G(T)  --PRECISION 15E.16
C15   G = ((((((((((((+3.620864000000000D-14)*T-1.051545600000000D-13)
C15  *    *T+2.245990400000000D-13)*T-8.785049600000000D-13)
C15  *    *T+3.893016320000000D-12)*T-1.855351552000000D-11)
C15  *    *T+1.039448166400000D-10)*T-7.151018001600000D-10)
C15  *    *T+6.420133522640000D-9)*T-8.291960820844000D-8)
C15  *    *T+1.821201851167060D-6)*T-9.621458822050362D-5)*T +
C15  *    4.677687402744898D-2
C
C     EXPANSION (0034) EVALUATED AS G(T)  --PRECISION 17E.18
      G = (((((((((((((((-2.29376000000000000D-15)
     *    *T+5.32480000000000000D-15)*T-4.83328000000000000D-15)
     *    *T+1.75718400000000000D-14)*T-7.43936000000000000D-14)
     *    *T+2.50224640000000000D-13)*T-9.23228160000000000D-13)
     *    *T+3.87554432000000000D-12)*T-1.85248000000000000D-11)
     *    *T+1.03950931840000000D-10)*T-7.15110504800000000D-10)
     *    *T+6.42013250344000000D-9)*T-8.29196070929200000D-8)
     *    *T+1.82120185123076000D-6)*T-9.62145882205441600D-5)*T +
     *    4.67768740274489776D-2
C
C     EXPANSION (0034) EVALUATED AS G(T)  --PRECISION 19E.20
C19   G = ((((((((((((((((-6.5536000000000000000D-16
C19  *    *T+9.8304000000000000000D-16)*T+4.9152000000000000000D-16)
C19  *    *T+1.3926400000000000000D-15)*T-9.7075200000000000000D-15)
C19  *    *T+2.3961600000000000000D-14)*T-6.9867520000000000000D-14)
C19  *    *T+2.4481792000000000000D-13)*T-9.2562176000000000000D-13)
C19  *    *T+3.8780787200000000000D-12)*T-1.8524081920000000000D-11)
C19  *    *T+1.0395028672000000000D-10)*T-7.1511061904000000000D-10)
C19  *    *T+6.4201325840800000000D-9)*T-8.2919607084760000000D-8)
C19  *    *T+1.8212018512269200000D-6)*T-9.6214588220544330000D-5)*T +
C19  *    4.6776874027448977670D-2
C
      S17AFF = Y*CX + G*SX
      GO TO 80
C     ERROR 1 EXIT
   60 IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      S17AFF = SQRT(TBPI/T)
C
   80 RETURN
      END
