      SUBROUTINE D02QFU(IREVCM,K,GWANT,NEQG,KROOT,INROOT,GOLD,PROOT,
     *                  ROOTD,GP,IGSC,TOUT)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C
C     DATE WRITTEN   840908   (YYMMDD)
C     REVISION DATE  850101   (YYMMDD)
C     AUTHOR  WATTS, H. A., (SNLA)
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  GWANT, TOUT
      INTEGER           INROOT, IREVCM, K, KROOT, NEQG
C     .. Array Arguments ..
      DOUBLE PRECISION  GOLD(NEQG), GP(NEQG), PROOT(NEQG), ROOTD(NEQG)
      INTEGER           IGSC(NEQG)
C     .. Scalars in Common ..
      DOUBLE PRECISION  DELSGN, FOURU, SRBIG, SRU, SVTNEW, TLEFT,
     *                  TROOTS, U, U34, U78, ZERO
      INTEGER           INROTP, KROO, KROOTP
      LOGICAL           DISCOP, GSTOP, NEEDG, NEWGEQ, PGSTOP, PSERCH,
     *                  ROOT, ROOTS, SEARCH
C     .. Local Scalars ..
      DOUBLE PRECISION  GTST
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN
C     .. Common blocks ..
      COMMON            /BD02QF/ZERO, U, FOURU, SRU, U34, U78, SRBIG,
     *                  DELSGN, TROOTS, TLEFT, SVTNEW, KROOTP, INROTP,
     *                  GSTOP, PGSTOP, ROOT, ROOTS, NEEDG, DISCOP,
     *                  NEWGEQ, SEARCH, PSERCH
      COMMON            /HD02QF/KROO
C     .. Save statement ..
      SAVE              /BD02QF/, /HD02QF/
C     .. Executable Statements ..
C
C     -- IT IS PRESUMED THAT THE USER HAS TOLD THE CODE THAT INTEGRATION
C        CANNOT BE CARRIED OUT BEYOND SOME POINT TSTOP. SECONDLY, IN
C        STEPPING TO TSTOP, IT IS PRESUMED THAT THE INTERVAL IS TOO
C        SMALL TO USE THE BASIC INTEGRATION METHOD OVER. INSTEAD, AN
C        EXTRAPOLATION OF THE SOLUTION IS ASSUMED BY SOME MEANS, SUCH AS
C        USING THE EULER FORMULA. IN THESE SPECIAL CIRCUMSTANCES, THIS
C        ROUTINE CAN BE CALLED TO CHECK FOR A ROOT ON THIS EXTREMELY
C        TINY INTERVAL.
C
      IF (IREVCM.EQ.12) GO TO 20
      K = 1
      IREVCM = 12
      RETURN
   20 CONTINUE
C     DO 20 K=1,NEQG
C       GTST = GRF(TOUT,Y,YPOUT,K,RPAR,IPAR)
      GTST = GWANT
      IF (ABS(GTST).GT.ZERO) GO TO 40
      IF (ABS(IGSC(K)).NE.1) IGSC(K) = 2
      GO TO 60
   40 IF (SIGN(1.D0,GOLD(K)).EQ.SIGN(1.D0,GTST)) GO TO 100
      IGSC(K) = 1
      IF (GTST.LT.0.D0) IGSC(K) = -1
   60 ROOT = .TRUE.
      KROO = K
      GOLD(K) = GTST
      GP(K) = GTST
      IF ( .NOT. SEARCH) GO TO 80
      IF (PROOT(K).NE.SRBIG) ROOTD(K) = ABS(PROOT(K)-TOUT)
   80 PROOT(K) = TOUT
C     20 CONTINUE
  100 K = K + 1
      IF (K.LE.NEQG) RETURN
      IREVCM = 0
      IF ( .NOT. ROOT) RETURN
      ROOTS = .TRUE.
      TROOTS = TOUT
      KROOTP = KROOT
      KROOT = KROO
      INROTP = INROOT
      INROOT = 1
      RETURN
C
C
C     END OF D02QFU (RDEO)
C
C
      END
