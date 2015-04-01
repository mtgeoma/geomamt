      INTEGER FUNCTION G05DRF(ALAMDA,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     G05DRF returns a Poisson variate from a distribution with the
C     supplied parameter. It returns zero if the supplied parameter
C     is negative or greater than MAXINT/2.
C
C
C
C
C
C
C
C     .. Parameters ..
      CHARACTER*6             SRNAME
      PARAMETER               (SRNAME='G05DRF')
      DOUBLE PRECISION        INVSQR
      PARAMETER               (INVSQR=0.39894228040143267794D0)
      DOUBLE PRECISION        ZERO, QUART, HALF, ONE, TWO
      PARAMETER               (ZERO=0.0D0,QUART=0.25D0,HALF=0.5D0,
     *                        ONE=1.0D0,TWO=2.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION        ALAMDA
      INTEGER                 IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION        D, INVARG, INVEXP, K, LASTAR, P, Q, R, S,
     *                        SQINAR, T, T1, T2, T3, U, V, W1, WR, WZ,
     *                        X, Y, Z, ZEROPR
      INTEGER                 I, IERROR, N, NREC
C     .. Local Arrays ..
      DOUBLE PRECISION        FAC(171)
      CHARACTER*80            P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION        G05CAF
      INTEGER                 P01ABF, X02BBF
      EXTERNAL                G05CAF, P01ABF, X02BBF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, EXP, INT, LOG, MOD, DBLE, SQRT
C     .. Save statement ..
      SAVE                    LASTAR, INVARG, INVEXP, SQINAR, ZEROPR, U,
     *                        V, T1, T2, T3, W1, WR, WZ
C     .. Data statements ..
C
      DATA                    LASTAR/-1.0D0/
      DATA                    FAC(1)/1.00000000000000000000D0/,
     *                        FAC(2)/0.70710678118654752440D0/,
     *                        FAC(3)/0.55032120814910444731D0/,
     *                        FAC(4)/0.45180100180492241598D0/,
     *                        FAC(5)/0.38385194963737748770D0/,
     *                        FAC(6)/0.33402418826640123076D0/,
     *                        FAC(7)/0.29585666126244628008D0/,
     *                        FAC(8)/0.26565006993025408476D0/,
     *                        FAC(9)/0.24112850410016933670D0/,
     *                        FAC(10)/0.22081252132060088614D0/
      DATA                    FAC(11)/0.20369756797298371676D0/,
     *                        FAC(12)/0.18907694923493627077D0/,
     *                        FAC(13)/0.17643868885977913751D0/,
     *                        FAC(14)/0.16540257277869528339D0/,
     *                        FAC(15)/0.15568019225792079758D0/,
     *                        FAC(16)/0.14704872910411693759D0/,
     *                        FAC(17)/0.13933325831461853748D0/,
     *                        FAC(18)/0.13239449985483095327D0/,
     *                        FAC(19)/0.12612015436168651333D0/,
     *                        FAC(20)/0.12041865418254125456D0/
      DATA                    FAC(21)/0.11521457781819726253D0/,
     *                        FAC(22)/0.11044523232139428893D0/,
     *                        FAC(23)/0.10605807018446566882D0/,
     *                        FAC(24)/0.10200871193123156262D0/,
     *                        FAC(25)/0.09825941469860363433D0/,
     *                        FAC(26)/0.09477787353433615834D0/,
     *                        FAC(27)/0.09153627390427658080D0/,
     *                        FAC(28)/0.08851053597663465557D0/,
     *                        FAC(29)/0.08567970681319882297D0/,
     *                        FAC(30)/0.08302546771642704789D0/
      DATA                    FAC(31)/0.08053173202435457627D0/,
     *                        FAC(32)/0.07818431453033058356D0/,
     *                        FAC(33)/0.07597065805690013658D0/,
     *                        FAC(34)/0.07387960596411258998D0/,
     *                        FAC(35)/0.07190121182358203980D0/,
     *                        FAC(36)/0.07002657935378538321D0/,
     *                        FAC(37)/0.06824772714156463943D0/,
     *                        FAC(38)/0.06655747377945760923D0/,
     *                        FAC(39)/0.06494933990836415693D0/,
     *                        FAC(40)/0.06341746432902251740D0/
      DATA                    FAC(41)/0.06195653187746719162D0/,
     *                        FAC(42)/0.06056171118169409288D0/,
     *                        FAC(43)/0.05922860075374567552D0/,
     *                        FAC(44)/0.05795318214199474137D0/,
     *                        FAC(45)/0.05673177908679455700D0/,
     *                        FAC(46)/0.05556102179982744928D0/,
     *                        FAC(47)/0.05443781563189839402D0/,
     *                        FAC(48)/0.05335931351217985431D0/,
     *                        FAC(49)/0.05232289163918178468D0/,
     *                        FAC(50)/0.05132612798405901291D0/
      DATA                    FAC(51)/0.05036678323349091964D0/,
     *                        FAC(52)/0.04944278385483029306D0/,
     *                        FAC(53)/0.04855220701256009781D0/,
     *                        FAC(54)/0.04769326710395605104D0/,
     *                        FAC(55)/0.04686430371454890626D0/,
     *                        FAC(56)/0.04606377082158110511D0/,
     *                        FAC(57)/0.04529022709702569483D0/,
     *                        FAC(58)/0.04454232718158894009D0/,
     *                        FAC(59)/0.04381881381803142421D0/,
     *                        FAC(60)/0.04311851074659185835D0/
      DATA                    FAC(61)/0.04244031627767598768D0/,
     *                        FAC(62)/0.04178319746760479337D0/,
     *                        FAC(63)/0.04114618483237131441D0/,
     *                        FAC(64)/0.04052836754225839685D0/,
     *                        FAC(65)/0.03992888904700805046D0/,
     *                        FAC(66)/0.03934694308716393254D0/,
     *                        FAC(67)/0.03878177005236372994D0/,
     *                        FAC(68)/0.03823265365184900826D0/,
     *                        FAC(69)/0.03769891786638030172D0/,
     *                        FAC(70)/0.03717992415417429773D0/
      DATA                    FAC(71)/0.03667506888648538710D0/,
     *                        FAC(72)/0.03618378099109300148D0/,
     *                        FAC(73)/0.03570551978427800225D0/,
     *                        FAC(74)/0.03523977297391776799D0/,
     *                        FAC(75)/0.03478605481813640733D0/,
     *                        FAC(76)/0.03434390442554448468D0/,
     *                        FAC(77)/0.03391288418451827666D0/,
     *                        FAC(78)/0.03349257831022471270D0/,
     *                        FAC(79)/0.03308259149921452359D0/,
     *                        FAC(80)/0.03268254768239980146D0/
      DATA                    FAC(81)/0.03229208886811797168D0/,
     *                        FAC(82)/0.03191087406777495668D0/,
     *                        FAC(83)/0.03153857829726726842D0/,
     *                        FAC(84)/0.03117489164801566386D0/,
     *                        FAC(85)/0.03081951842201035886D0/,
     *                        FAC(86)/0.03047217632577707991D0/,
     *                        FAC(87)/0.03013259571863097280D0/,
     *                        FAC(88)/0.02980051891099731727D0/,
     *                        FAC(89)/0.02947569950894915734D0/,
     *                        FAC(90)/0.02915790180144678387D0/
      DATA                    FAC(91)/0.02884690018706641077D0/,
     *                        FAC(92)/0.02854247863727882119D0/,
     *                        FAC(93)/0.02824443019358627757D0/,
     *                        FAC(94)/0.02795255649605029162D0/,
     *                        FAC(95)/0.02766666734094633265D0/,
     *                        FAC(96)/0.02738658026546634092D0/,
     *                        FAC(97)/0.02711212015755789715D0/,
     *                        FAC(98)/0.02684311888914176066D0/,
     *                        FAC(99)/0.02657941497108872601D0/,
     *                        FAC(100)/0.02632085322846369736D0/
      DATA                    FAC(101)/0.02606728449466073480D0/,
     *                        FAC(102)/0.02581856532315865524D0/,
     *                        FAC(103)/0.02557455771572352693D0/,
     *                        FAC(104)/0.02533512886597293818D0/,
     *                        FAC(105)/0.02510015091729801447D0/,
     *                        FAC(106)/0.02486950073421349502D0/,
     *                        FAC(107)/0.02464305968627437862D0/,
     *                        FAC(108)/0.02442071344376026842D0/,
     *                        FAC(109)/0.02420235178438608899D0/,
     *                        FAC(110)/0.02398786841035076847D0/
      DATA                    FAC(111)/0.02377716077508418551D0/,
     *                        FAC(112)/0.02357012991909754177D0/,
     *                        FAC(113)/0.02336668031438367346D0/,
     *                        FAC(114)/0.02316671971685195915D0/,
     *                        FAC(115)/0.02297015902631769290D0/,
     *                        FAC(116)/0.02277691215359831826D0/,
     *                        FAC(117)/0.02258689589429898610D0/,
     *                        FAC(118)/0.02240002980889771298D0/,
     *                        FAC(119)/0.02221623610876616209D0/,
     *                        FAC(120)/0.02203543954778591621D0/
      DATA                    FAC(121)/0.02185756731924221496D0/,
     *                        FAC(122)/0.02168254895769762863D0/,
     *                        FAC(123)/0.02151031624556716508D0/,
     *                        FAC(124)/0.02134080312413397270D0/,
     *                        FAC(125)/0.02117394560876121748D0/,
     *                        FAC(126)/0.02100968170807097379D0/,
     *                        FAC(127)/0.02084795134687516619D0/,
     *                        FAC(128)/0.02068869629265681449D0/,
     *                        FAC(129)/0.02053186008541214259D0/,
     *                        FAC(130)/0.02037738797067558050D0/
      DATA                    FAC(131)/0.02022522683556038306D0/,
     *                        FAC(132)/0.02007532514765756461D0/,
     *                        FAC(133)/0.01992763289664516037D0/,
     *                        FAC(134)/0.01978210153846852105D0/,
     *                        FAC(135)/0.01963868394196047192D0/,
     *                        FAC(136)/0.01949733433777776317D0/,
     *                        FAC(137)/0.01935800826953734216D0/,
     *                        FAC(138)/0.01922066254704262617D0/,
     *                        FAC(139)/0.01908525520149617766D0/,
     *                        FAC(140)/0.01895174544260101384D0/
      DATA                    FAC(141)/0.01882009361745824543D0/,
     *                        FAC(142)/0.01869026117117386107D0/,
     *                        FAC(143)/0.01856221060909227850D0/,
     *                        FAC(144)/0.01843590546057879165D0/,
     *                        FAC(145)/0.01831131024427727550D0/,
     *                        FAC(146)/0.01818839043477348619D0/,
     *                        FAC(147)/0.01806711243059802901D0/,
     *                        FAC(148)/0.01794744352350657899D0/,
     *                        FAC(149)/0.01782935186897824103D0/,
     *                        FAC(150)/0.01771280645787604367D0/
      DATA                    FAC(151)/0.01759777708921648553D0/,
     *                        FAC(152)/0.01748423434399780682D0/,
     *                        FAC(153)/0.01737214956003925274D0/,
     *                        FAC(154)/0.01726149480778604016D0/,
     *                        FAC(155)/0.01715224286703704332D0/,
     *                        FAC(156)/0.01704436720455438820D0/,
     *                        FAC(157)/0.01693784195251619524D0/,
     *                        FAC(158)/0.01683264188777564593D0/,
     *                        FAC(159)/0.01672874241189137580D0/,
     *                        FAC(160)/0.01662611953189592220D0/
      DATA                    FAC(161)/0.01652474984177058615D0/,
     *                        FAC(162)/0.01642461050459660856D0/,
     *                        FAC(163)/0.01632567923535401867D0/,
     *                        FAC(164)/0.01622793428434089088D0/,
     *                        FAC(165)/0.01613135442118705020D0/,
     *                        FAC(166)/0.01603591891943750124D0/,
     *                        FAC(167)/0.01594160754168202413D0/,
     *                        FAC(168)/0.01584840052520848778D0/,
     *                        FAC(169)/0.01575627856815847986D0/,
     *                        FAC(170)/0.01566522281616484641D0/,
     *                        FAC(171)/0.01557521484945167664D0/
C     .. Executable Statements ..
C
      IERROR = 0
      NREC = 0
C
C     Deal with invalid input parameter
C
      IF (ALAMDA.LE.ZERO) THEN
         IERROR = 1
         NREC = 1
         G05DRF = 0
         WRITE (P01REC,FMT=99999) ALAMDA
      ELSE IF (ALAMDA.GT.DBLE(X02BBF(0.0D0))/2) THEN
         IERROR = 2
         NREC = 1
         G05DRF = 0
         WRITE (P01REC,FMT=99998) ALAMDA
C
C     Deal with small input parameter by the product of uniforms method
C
      ELSE IF (ALAMDA.LT.7.50D0) THEN
         IF (ALAMDA.NE.LASTAR) THEN
            LASTAR = ALAMDA
            INVEXP = EXP(-ALAMDA)
            Z = INVEXP
         ELSE
            Z = INVEXP
         END IF
         N = 0
         X = ONE
   20    X = X*G05CAF(X)
         IF (X.GE.Z) THEN
            N = N + 1
            GO TO 20
         END IF
         G05DRF = N
      ELSE
C
C     Set up the contants for the acceptance and rejection functions
C
         IF (ALAMDA.NE.LASTAR) THEN
            LASTAR = ALAMDA
            INVARG = ONE/ALAMDA
            SQINAR = SQRT(INVARG)
            ZEROPR = ZERO
            U = ALAMDA + HALF
            V = 1.051D0/SQINAR + 0.305D0 + 0.2D0*SQINAR
            T1 = ONE + 0.0435D0*INVARG
            T2 = T1*(ONE-0.13D0*INVARG)
            T3 = T1*(ONE+QUART*SQINAR)
            W1 = 0.475D0 + 0.413D0*SQINAR + 0.75D0*INVARG
            WR = 0.475D0 + 0.162D0*SQINAR + 0.35D0*INVARG
            WZ = 0.552D0
         END IF
C
C     Go round in a loop generating new numbers and rejecting
C     ridiculous cases
C
   40    CONTINUE
         X = TWO*G05CAF(X) - ONE
         Y = TWO*G05CAF(Y) - ONE
         IF (X.EQ.ZERO .AND. Y.EQ.ZERO) GO TO 40
         IF (ABS(X).LE.ABS(Y)) THEN
            X = X/(Y*Y)
         ELSE
            X = Y/(X*X)
         END IF
         Y = U + V*X
         IF (Y.LT.ZERO) GO TO 40
         N = INT(Y)
         Y = X + HALF/V
C
C     Perform the quick acceptance test
C
         IF (Y.LE.ZERO) THEN
            P = ONE + W1*Y
         ELSE
            P = ONE - WR*Y
         END IF
         P = T2*P
         IF (ABS(X).LE.ONE) THEN
            Q = T1
         ELSE
            Q = T1/ABS(X*X*X)
         END IF
         R = Q*G05CAF(R)
         IF (R.LE.P) THEN
            G05DRF = N
            GO TO 120
         END IF
         Y = WZ*Y*Y
         P = T3 - Y + HALF*Y*Y
         IF (R.GT.P) GO TO 40
C
C     If the input parameter is fairly small, calculate the probability
C     directly. The complication is for speed and to avoid overflow.
C     Note that the probability of the result exceeding 171 is less than
C     1e-15 even when the parameter is 87, so the truncation is reasonab
C
         IF (ALAMDA.LT.87.0D0) THEN
            IF (ZEROPR.NE.ZERO) THEN
               P = ZEROPR
            ELSE
               ZEROPR = (-ALAMDA) - LOG(INVSQR*SQINAR)
               P = ZEROPR
            END IF
            IF (N.GT.0) THEN
               IF (N.LT.171) THEN
                  X = LOG(ALAMDA*FAC(N))
                  I = N
   60             CONTINUE
                  IF (MOD(I,2).NE.0) P = P + X
                  I = I/2
                  IF (I.EQ.0) GO TO 100
                  X = X + X
                  GO TO 60
               ELSE
                  GO TO 40
               END IF
            END IF
C
C     If the parameter is large, calculate the probability to 1e-15
C     accuracy (if the arithmetic permits) using Stirling's formula.
C     Note that the probability of the result being below 24 is less
C     than 1e-15 even when the parameter is 87. The complexity is to
C     keep the cancellation errors down to max(4096,sqrt(alamda))
C     times the machine precision.
C
         ELSE
            IF (N.LE.0) GO TO 40
            X = DBLE(N)
            D = (X-ALAMDA)/ALAMDA
            IF (ABS(D).LE.0.022D0) THEN
               K = ONE
               S = -(X-ALAMDA)*D
               T = -HALF*D
               Z = ZERO
   80          CONTINUE
               Y = S/(K*(K+ONE))
               Z = Z + Y + T/K
               K = K + ONE
               S = -D*S
               T = -D*T
               IF (ABS(Y).GT.1.0D-15) GO TO 80
            ELSE
               Z = (X-ALAMDA) - (X+HALF)*LOG(X/ALAMDA)
            END IF
            Y = ONE/(X*X)
            P = (Z-(0.08333333333333333333D0-Y*
     *          (0.002777777777777777778D0-Y*
     *          (0.0007936507936507936508D0-Y*
     *          0.0005952380952380952380D0)))/X)
         END IF
  100    IF (LOG(R).LE.P) THEN
            G05DRF = N
         ELSE
            GO TO 40
         END IF
      END IF
  120 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (1X,'** On entry, ALAMDA.le.ZERO : ALAMDA = ',D13.5)
99998 FORMAT (1X,'** On entry, 2*ALAMDA.gt.MAXINT : ALAMDA = ',D13.5)
      END
