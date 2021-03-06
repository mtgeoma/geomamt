      DOUBLE PRECISION FUNCTION S10AAF(X,IFAIL)
C     MARK 5A REVISED - NAG COPYRIGHT 1976
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     TANH(X)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 T, XHI, Y
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, EXP, SIGN
C     .. Data statements ..
C     PRECISION DEPENDENT CONSTANTS
C08   DATA XHI/9.25D0/
C12   DATA XHI/13.875D0/
C14   DATA XHI/16.25D0/
      DATA XHI/18.5D0/
C18   DATA XHI/20.75D0/
C     .. Executable Statements ..
C
C     NO FAILURE EXITS
      IFAIL = 0
C
C     TEST FOR ARGUMENT IN SMALLEST RANGE
      IF (ABS(X).GT.1.0D0) GO TO 20
C
C     ARGUMENT OF EXPANSION
      T = 2.0D0*X*X - 1.0D0
C
C      * EXPANSION (0010) *
C
C     EXPANSION (0010) EVALUATED AS Y(T)  --PRECISION 08E
C08   Y = (((((((-2.7725267D-6)*T+1.6336774D-5)*T-9.1410867D-5)
C08  *    *T+5.4272941D-4)*T-3.2252921D-3)*T+1.9180840D-2)
C08  *    *T-1.1588345D-1)*T + 8.6105717D-1
C
C     EXPANSION (0010) EVALUATED AS Y(T)  --PRECISION 12E
C12   Y = (((((((((((-2.29992850509D-9)*T+1.35520400383D-8)
C12  *    *T-7.35288916190D-8)*T+4.36647778988D-7)
C12  *    *T-2.59918065059D-6)*T+1.54253634350D-5)
C12  *    *T-9.15428534421D-5)*T+5.43306980751D-4)
C12  *    *T-3.22525518655D-3)*T+1.91807238533D-2)
C12  *    *T-1.15883448973D-1)*T + 8.61057171580D-1
C
C     EXPANSION (0010) EVALUATED AS Y(T)  --PRECISION 14E
C14   Y = (((((((((((((-6.6242033444587D-11)*T+3.9032286757454D-10)
C14  *    *T-2.0846418963935D-9)*T+1.2381071435532D-8)
C14  *    *T-7.3797999879855D-8)*T+4.3796511866652D-7)
C14  *    *T-2.5990191856299D-6)*T+1.5424680370018D-5)
C14  *    *T-9.1542900536031D-5)*T+5.4330714084406D-4)
C14  *    *T-3.2252551806632D-3)*T+1.9180723839605D-2)
C14  *    *T-1.1588344897285D-1)*T + 8.6105717158055D-1
C
C     EXPANSION (0010) EVALUATED AS Y(T)  --PRECISION 16E
      Y = (((((((((((((((-1.907888434471600D-12)
     *    *T+1.124199312776748D-11)*T-5.908745181531817D-11)
     *    *T+3.509758916273561D-10)*T-2.095373768837420D-9)
     *    *T+1.243517352745986D-8)*T-7.378980192173815D-8)
     *    *T+4.379282308765732D-7)*T-2.599022539340038D-6)
     *    *T+1.542469328074432D-5)*T-9.154289983175165D-5)
     *    *T+5.433071386922689D-4)*T-3.225255180728459D-3)
     *    *T+1.918072383973950D-2)*T-1.158834489728470D-1)*T +
     *     8.610571715805476D-1
C
C     EXPANSION (0010) EVALUATED AS Y(T)  --PRECISION 18E
C18   Y = ((((((((((((((((+3.23789406115168256D-13)
C18  *    *T-1.90788843447160013D-12)*T+9.94683550330680934D-12)
C18  *    *T-5.90874518153181676D-11)*T+3.53080522767104662D-10)
C18  *    *T-2.09537376883742008D-9)*T+1.24333926857262297D-8)
C18  *    *T-7.37898019217381456D-8)*T+4.37929065646135880D-7)
C18  *    *T-2.59902253934003830D-6)*T+1.54246930682575238D-5)
C18  *    *T-9.15428998317516507D-5)*T+5.43307138718829741D-4)
C18  *    *T-3.22525518072845854D-3)*T+1.91807238397382333D-2)
C18  *    *T-1.15883448972846974D-1)*T + 8.61057171580547644D-1
C
C
      S10AAF = X*Y
      RETURN
C
C     TEST FOR ARGUMENT IN MIDDLE RANGE
   20 IF (ABS(X).GE.XHI) GO TO 40
      Y = EXP(2.0D0*X)
      S10AAF = (Y-1.0D0)/(Y+1.0D0)
      RETURN
C
C     LARGE ARGUMENTS
   40 S10AAF = SIGN(1.0D0,X)
      RETURN
C
      END
