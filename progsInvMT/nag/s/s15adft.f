      DOUBLE PRECISION FUNCTION S15ADF(X,IFAIL)
C     MARK 5A REVISED - NAG COPYRIGHT 1976
C     MARK 5C REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     COMPLEMENT OF ERROR FUNCTION ERFC(X)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 T, XHI, XLO, Y
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, EXP
C     .. Data statements ..
C     PRECISION DEPENDENT CONSTANTS
C08   DATA XLO/-4.5D0/
C12   DATA XLO/-5.25D0/
C14   DATA XLO/-5.75D0/
      DATA XLO/-6.25D0/
C18   DATA XLO/-6.5D0/
C
C     RANGE DEPENDENT CONSTANTS
      DATA XHI/ 2.66D+1 /
C     XHI = LARGEST X SUCH THAT EXP(-X*X) .GT. MINREAL (ROUNDED DOWN)
CR1   DATA XHI/13.0D0/
CR2   DATA XHI/9.5D0/
CR3   DATA XHI/13.0D0/
CR4   DATA XHI/25.0D0/
CR5   DATA XHI/26.0D0/
C     .. Executable Statements ..
C
C     NO FAILURE EXITS
      IFAIL = 0
C     TEST EXTREME EXITS
      IF (X.GE.XHI) GO TO 20
      IF (X.LE.XLO) GO TO 40
C
C     EXPANSION ARGUMENT
      T = 1.0D0 - 7.5D0/(ABS(X)+3.75D0)
C
C      * EXPANSION (0021) *
C
C     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 08E
C08   Y = (((((((((((+3.1475326D-5)*T-1.3874589D-4)*T-6.4127909D-6)
C08  *    *T+1.7866301D-3)*T-8.2316935D-3)*T+2.4151896D-2)
C08  *    *T-5.4799165D-2)*T+1.0260225D-1)*T-1.6357229D-1)
C08  *    *T+2.2600824D-1)*T-2.7342192D-1)*T + 1.4558972D-1
C
C     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 12E
C12   Y = ((((((((((((((((-4.21661579602D-8*T-8.63384346353D-8)
C12  *    *T+6.06038693567D-7)*T+5.90655413508D-7)
C12  *    *T-6.12872971594D-6)*T+3.73223486059D-6)
C12  *    *T+4.78645837248D-5)*T-1.52546487034D-4)
C12  *    *T-2.55222360474D-5)*T+1.80299061562D-3)
C12  *    *T-8.22062412199D-3)*T+2.41432185990D-2)
C12  *    *T-5.48023263289D-2)*T+1.02604312548D-1)
C12  *    *T-1.63571895545D-1)*T+2.26008066898D-1)
C12  *    *T-2.73421931495D-1)*T + 1.45589721275D-1
C
C     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 14E
C14   Y = (((((((((((((((-2.2356173494379D-9
C14  *    *T+4.5302502889845D-9)*T+2.5918103316137D-8)
C14  *    *T-6.3684846832869D-8)*T-1.7642194353331D-7)
C14  *    *T+6.4907607131235D-7)*T+7.4296952017617D-7)
C14  *    *T-6.1758018478516D-6)*T+3.5866167916231D-6)
C14  *    *T+4.7895180610590D-5)*T-1.5246364229106D-4)
C14  *    *T-2.5534256252531D-5)*T+1.8029626230333D-3)
C14  *    *T-8.2206213481002D-3)*T+2.4143223946968D-2)
C14  *    *T-5.4802326675661D-2)*T+1.0260431203382D-1
C14   Y = (((Y*T-1.6357189552481D-1)*T+2.2600806691658D-1)
C14  *    *T-2.7342193149541D-1)*T + 1.4558972127504D-1
C
C     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 16E
      Y = (((((((((((((((+3.328130055126039D-10
     *    *T-5.718639670776992D-10)*T-4.066088879757269D-9)
     *    *T+7.532536116142436D-9)*T+3.026547320064576D-8)
     *    *T-7.043998994397452D-8)*T-1.822565715362025D-7)
     *    *T+6.575825478226343D-7)*T+7.478317101785790D-7)
     *    *T-6.182369348098529D-6)*T+3.584014089915968D-6)
     *    *T+4.789838226695987D-5)*T-1.524627476123466D-4)
     *    *T-2.553523453642242D-5)*T+1.802962431316418D-3)
     *    *T-8.220621168415435D-3)*T+2.414322397093253D-2
      Y = (((((Y*T-5.480232669380236D-2)*T+1.026043120322792D-1)
     *    *T-1.635718955239687D-1)*T+2.260080669166197D-1)
     *    *T-2.734219314954260D-1)*T + 1.455897212750385D-1
C
C     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 18E
C18   Y = (((((((((((((((-1.58023488119651697D-11
C18  *    *T-4.94972069009392927D-11)*T+1.86424953544623784D-10)
C18  *    *T+6.29796246918239617D-10)*T-1.34751340973493898D-9)
C18  *    *T-4.84566988844706300D-9)*T+9.22474802259858004D-9)
C18  *    *T+3.14410318645430670D-8)*T-7.26754673242913196D-8)
C18  *    *T-1.83380699508554268D-7)*T+6.59488268069175234D-7)
C18  *    *T+7.48541685740064308D-7)*T-6.18344429012694168D-6)
C18  *    *T+3.58371497984145357D-6)*T+4.78987832434182054D-5)
C18  *    *T-1.52462664665855354D-4)*T-2.55353311432760448D-5
C18   Y = ((((((((Y*T+1.80296241673597993D-3)
C18  *    *T-8.22062115413991215D-3)
C18  *    *T+2.41432239724445769D-2)*T-5.48023266949776152D-2)
C18  *    *T+1.02604312032198239D-1)*T-1.63571895523923969D-1)
C18  *    *T+2.26008066916621431D-1)*T-2.73421931495426482D-1)*T +
C18  *     1.45589721275038539D-1
C
      S15ADF = EXP(-X*X)*Y
      IF (X.LT.0.0D0) S15ADF = 2.0D0 - S15ADF
      RETURN
C
   20 S15ADF = 0.0D0
      RETURN
   40 S15ADF = 2.0D0
      RETURN
C
      END
