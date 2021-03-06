      DOUBLE PRECISION FUNCTION S14BAX(X)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     Computes the value X - log(1+X), retaining full relative
C     precision. Uses two Chebyshev expansions, the first valid in the
C     range -0.5 .lt. X .lt. 0.0, the second valid in the range
C     0.0 .lt. X .lt. 0.5.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. Local Scalars ..
      DOUBLE PRECISION                 Z
C     .. Executable Statements ..
      IF (X.GT.0.0D0) THEN
         Z = 4*X - 1
C        Precision           20 sig. figs.
         S14BAX = ((((((((((((((-2.967444749673129063172D-15+
     *            (5.727889008167279942775D-16)*Z)
     *            *Z+1.254316797763179284544D-14)
     *            *Z+(-6.609208145248494083774D-14))
     *            *Z+3.552093631141621449360D-13)
     *            *Z+(-1.877575086584465734489D-12))
     *            *Z+9.950448268789644316317D-12)
     *            *Z+(-5.297680484461551503860D-11))
     *            *Z+2.832511945923442439052D-10)
     *            *Z+(-1.521770786533239703319D-09))
     *            *Z+8.221622611294801319744D-09)
     *            *Z+(-4.471056587610118468800D-08))
     *            *Z+2.450395094674988857885D-07)
     *            *Z+(-1.355590675105630409968D-06))
     *            *Z+7.586141840685986975605D-06)*Z
         S14BAX = X*X*((((((S14BAX+(-4.307383586344297282484D-05))
     *            *Z+2.492281965528721274994D-04)
     *            *Z+(-1.479382557242297784497D-03))
     *            *Z+9.109536917931723215386D-03)
     *            *Z+(-5.940635794528781547636D-02))
     *            *Z+4.297031789726439077393D-01)
      ELSE
         Z = 4*X + 1
C        Precision           20 sig. figs.
         S14BAX = ((((((((((((((-9.222006652207022772023D-14+
     *            (3.049209823713897715108D-14)*Z)
     *            *Z+8.110971610764310768061D-14)
     *            *Z+(-2.708761967476698009188D-13))
     *            *Z+1.468315445990726325182D-12)
     *            *Z+(-4.549612295389177687275D-12))
     *            *Z+1.318723069570494295267D-11)
     *            *Z+(-4.157371211079606768770D-11))
     *            *Z+1.323612522929628013599D-10)
     *            *Z+(-4.184113387833978015279D-10))
     *            *Z+1.325817212562882710866D-09)
     *            *Z+(-4.217964594365521564431D-09))
     *            *Z+1.346846941394504845741D-08)
     *            *Z+(-4.318289149822655693289D-08))
     *            *Z+1.391089264405774451800D-07)*Z
         S14BAX = X*X*((((((((((((S14BAX+(-4.505690978103018400967D-07))
     *            *Z+1.468654737538004800256D-06)
     *            *Z+(-4.823073082193340091105D-06))
     *            *Z+1.598133959833499951509D-05)
     *            *Z+(-5.353471603473257684307D-05))
     *            *Z+1.817808088831846419015D-04)
     *            *Z+(-6.280405138046454551685D-04))
     *            *Z+2.220117130128515512154D-03)
     *            *Z+(-8.100449505773729833165D-03))
     *            *Z+3.096169990770673931421D-02)
     *            *Z+(-1.275070148763436552823D-01))
     *            *Z+6.029131592284948390275D-01)
      END IF
      RETURN
      END
