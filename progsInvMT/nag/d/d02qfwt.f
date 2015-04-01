      SUBROUTINE D02QFW(NEQG,KROOT,INROOT,TKT,GOLD,PROOT,ROOTD,GP,
     *                  NEEDGK,IGSC,T)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C
C     DATE WRITTEN   840908   (YYMMDD)
C     REVISION DATE  850101   (YYMMDD)
C     AUTHOR  WATTS, H. A., (SNLA)
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           INROOT, KROOT, NEQG
C     .. Array Arguments ..
      DOUBLE PRECISION  GOLD(NEQG), GP(NEQG), PROOT(NEQG), ROOTD(NEQG),
     *                  TKT(NEQG)
      INTEGER           IGSC(NEQG), NEEDGK(NEQG)
C     .. Scalars in Common ..
      DOUBLE PRECISION  DELSGN, FOURU, SRBIG, SRU, SVTNEW, TLEFT,
     *                  TROOTS, U, U34, U78, ZERO
      INTEGER           INROTP, KROOTP
      LOGICAL           DISCOP, GSTOP, NEEDG, NEWGEQ, PGSTOP, PSERCH,
     *                  ROOT, ROOTS, SEARCH
C     .. Local Scalars ..
      INTEGER           K
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /BD02QF/ZERO, U, FOURU, SRU, U34, U78, SRBIG,
     *                  DELSGN, TROOTS, TLEFT, SVTNEW, KROOTP, INROTP,
     *                  GSTOP, PGSTOP, ROOT, ROOTS, NEEDG, DISCOP,
     *                  NEWGEQ, SEARCH, PSERCH
C     .. Save statement ..
      SAVE              /BD02QF/
C     .. Executable Statements ..
C
C                       -- INITIALIZATION AT THE VERY BEGINNING OF THE
C                          ROOT SEARCHING PROCESS AND WHENEVER THE ROOT
C                          FUNCTIONS HAVE BEEN ALTERED
C                       -- EVALUATE GRF AND CHECK FOR A ROOT
C
      NEEDG = .TRUE.
      ROOT = .FALSE.
      ROOTS = .FALSE.
      KROOT = 0
      KROOTP = 0
      INROOT = 1
      INROTP = 0
      TLEFT = T
      SVTNEW = T
      DO 40 K = 1, NEQG
         TKT(K) = T
         IGSC(K) = 0
         PROOT(K) = SRBIG
         IF ( .NOT. SEARCH) GO TO 20
         ROOTD(K) = SRBIG
         NEEDGK(K) = 0
C        10   GOLD(K) = GRF(T,Y,YPOUT,K,RPAR,IPAR)
   20    GP(K) = GOLD(K)
         IF (ABS(GOLD(K)).GT.ZERO) GO TO 40
         ROOT = .TRUE.
         ROOTS = .TRUE.
         TROOTS = T
         KROOT = K
         IGSC(K) = 2
         PROOT(K) = T
   40 CONTINUE
      RETURN
C
C
C     END OF D02QFW (RDEI)
C
C
      END
