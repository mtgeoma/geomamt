      DOUBLE PRECISION FUNCTION S01BAF(X,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Computes log(1+x) retaining full precision when x is small.
C     Following a suggestion of N.M.Temme, uses the Chebyshev expansion
C       log((1+p*p+2*p*z)/(1+p*p-2*p*z)) =
C         4 * SUM (k=0 to infinity) p**(2*k+1) * T(2*k+1)(z) / (2*k+1)
C     (Lyusternik et al, Handbook for Computing Elementary Functions,
C     Pergamon Press, 1965).
C     Equating (1+x) with (1+p*p+2*p*z)/(1+p*p-2*p*z), we get
C       z = x*(1+p*p)/(2*p*(x+2)).
C     Choosing p as (sqrt(sqrt(2))-1)/(sqrt(sqrt(2))+1) =
C     0.0864..., the expansion will converge with a very small number
C     of terms, and be valid in the range [1/sqrt(2)-1,sqrt(2)-1].
C     Outside this range, the standard Fortran log function is good
C     enough.
C
C     **************************************************************
C
C     TO EXTRACT THE CORRECT CODE FOR A PARTICULAR MACHINE-RANGE,
C     ACTIVATE THE STATEMENTS CONTAINED IN COMMENTS BEGINNING  CDD ,
C     WHERE  DD  IS THE APPROXIMATE NUMBER OF SIGNIFICANT DECIMAL
C     DIGITS REPRESENTED BY THE MACHINE
C     DELETE THE ILLEGAL DUMMY STATEMENT OF THE FORM
C     * EXPANSION (NNNN) *
C
C     **************************************************************
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='S01BAF')
      DOUBLE PRECISION                 ZERO, ONE, TWO, FOUR
      PARAMETER                        (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,
     *                                 FOUR=4.0D0)
      DOUBLE PRECISION                 P, TOPCON, BOTCON
      PARAMETER                        (P=0.864272337258897920439D-1,
     *                                 TOPCON=ONE+P*P,BOTCON=TWO*P)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 XBAR, XBAR2, Y
      INTEGER                          IER, NREC
C     .. Local Arrays ..
      CHARACTER*80                     REC(1)
C     .. External Functions ..
      INTEGER                          P01ABF
      EXTERNAL                         P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        LOG
C     .. Executable Statements ..
      IER = 0
      NREC = 0
      IF (X.LE.-ONE) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT=99999) X
         S01BAF = ZERO
      ELSE
         XBAR = X*TOPCON/((X+TWO)*BOTCON)
         XBAR2 = XBAR*XBAR
         IF (XBAR2.GT.ONE) THEN
C           X lies outside the domain of the Chebyshev expansion. The
C           standard log function is good enough here.
            S01BAF = LOG(ONE+X)
         ELSE
C      * EXPANSION (1989) *
C
C           Precision 8 sig. figs.
C08         Y = ((((3.293337379D-07)*XBAR2+1.485498471D-05)
C08  *          *XBAR2+8.417758616D-04)*XBAR2+8.578643736D-02)*XBAR
C
C           Precision 9 sig. figs.
C09         Y = (((((7.65337460178D-09)*XBAR2+3.12113645041D-07)
C09  *          *XBAR2+1.48678997788D-05)*XBAR2+8.41772274111D-04)
C09  *          *XBAR2+8.57864376289D-02)*XBAR
C
C           Precision 12 sig. figs.
C12         Y = ((((((1.8709578861177D-10)*XBAR2+7.1388611830954D-09)
C12  *          *XBAR2+3.1262815845939D-07)*XBAR2+1.4867674679215D-05)
C12  *          *XBAR2+8.4177231430731D-04)*XBAR2+8.5786437626890D-02)
C12  *          *XBAR
C
C           Precision 15 sig. figs.
C15         Y = ((((((((1.22486400375197395D-13)
C15  *          *XBAR2+4.27082217141498675D-12)
C15  *          *XBAR2+1.72411799552206238D-10)
C15  *          *XBAR2+7.15755109317092088D-09)
C15  *          *XBAR2+3.12616844036222621D-07)
C15  *          *XBAR2+1.48676779968385945D-05)
C15  *          *XBAR2+8.41772313891137149D-04)
C15  *          *XBAR2+8.57864376269049504D-02)*XBAR
C
C           Precision 17 sig. figs.
            Y = (((((((((3.2291738460000068842D-15)
     *          *XBAR2+1.0876241152969736561D-13)
     *          *XBAR2+4.2948391518946117974D-12)
     *          *XBAR2+1.7238949807033230063D-10)
     *          *XBAR2+7.1575628872238349836D-09)
     *          *XBAR2+3.1261684049800674718D-07)
     *          *XBAR2+1.4867677997401492503D-05)
     *          *XBAR2+8.4177231389109694161D-04)
     *          *XBAR2+8.5786437626904951205D-02)*XBAR
C
C           Precision 19 sig. figs.
C19         Y = ((((((((((8.632726136875131478175D-17)
C19  *          *XBAR2+2.819119354498438138997D-15)
C19  *          *XBAR2+1.095825205127005030966D-13)
C19  *          *XBAR2+4.293942157694452115758D-12)
C19  *          *XBAR2+1.723900811165624044214D-10)
C19  *          *XBAR2+7.157562658169958871436D-09)
C19  *          *XBAR2+3.126168405508653339736D-07)
C19  *          *XBAR2+1.486767799739488518013D-05)
C19  *          *XBAR2+8.417723138910973170285D-04)
C19  *          *XBAR2+8.578643762690495119826D-02)*XBAR
C
            S01BAF = FOUR*Y
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, X.le.-1.0 : X = ',1P,D13.5)
      END
