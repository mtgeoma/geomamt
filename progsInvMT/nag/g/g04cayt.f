      SUBROUTINE G04CAY(NFAC,LFAC,NREPL,N,Y,R,E,TMEAN,IE,SEMEAN,IMEAN,
     *                  IR,IFAC,TABLE,LDT,ITERM,RSS,NRDF)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Sweep out main effects
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RSS
      INTEGER           IE, IR, ITERM, LDT, N, NFAC, NRDF, NREPL
C     .. Array Arguments ..
      DOUBLE PRECISION  E(*), R(N), SEMEAN(*), TABLE(LDT,5), TMEAN(*),
     *                  Y(N)
      INTEGER           IFAC(N), IMEAN(*), LFAC(NFAC)
C     .. Local Scalars ..
      DOUBLE PRECISION  ESS
      INTEGER           I, IEF, II, J, K, L, LFACI, NREPF, NREPI
C     .. External Subroutines ..
      EXTERNAL          G04CAW
C     .. Executable Statements ..
C
C     Initialise indicators
C
      ITERM = ITERM + NFAC
      IR = IR + NFAC
      IEF = IE
      DO 20 I = 1, NFAC
         IEF = IEF + LFAC(I)
   20 CONTINUE
      IE = IEF
C
C     NREPI = number of replicates for factor level
C     NREPF = number of times each complete set of factor levels is
C             repeated
C
      NREPI = NREPL
      DO 100 I = NFAC, 1, -1
         LFACI = LFAC(I)
         IMEAN(IR) = IE
         IE = IE - LFACI
         NREPF = N/(NREPI*LFACI)
         II = 0
         DO 80 J = 1, NREPF
            DO 60 K = 1, LFACI
               DO 40 L = 1, NREPI
                  II = II + 1
                  E(IE+K) = E(IE+K) + R(II)
                  TMEAN(IE+K) = TMEAN(IE+K) + Y(II)
                  IFAC(II) = K
   40          CONTINUE
   60       CONTINUE
   80    CONTINUE
         CALL G04CAW(N,LFACI,E(IE+1),TMEAN(IE+1),SEMEAN(IR),ESS,R,IFAC)
         TABLE(ITERM,1) = LFACI - 1
         TABLE(ITERM,2) = ESS
         NRDF = NRDF - LFACI + 1
         NREPI = NREPI*LFACI
         ITERM = ITERM - 1
         IR = IR - 1
  100 CONTINUE
C
C     Reset indicators
C
      ITERM = ITERM + NFAC
      IR = IR + NFAC
      IE = IEF
      RETURN
      END
