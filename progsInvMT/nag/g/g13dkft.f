      SUBROUTINE G13DKF(K,LMAX,M,MLAST,Z,IK,REF,LREF,V,PREDZ,SEFZ,WORK,
     *                  IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13DKF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IK, K, LMAX, LREF, M, MLAST
C     .. Array Arguments ..
      DOUBLE PRECISION  PREDZ(IK,LMAX), REF(LREF), SEFZ(IK,LMAX),
     *                  V(IK,M), WORK(K*M), Z(IK,M)
C     .. Local Scalars ..
      INTEGER           IERROR, JREF, LR1, LR2, LR3, LR4, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G13DKZ
C     .. Executable Statements ..
C
C     This subroutine updates a set of forecasts and their standard
C     deviations.  M additional observations are assumed available
C     for each of the k series. The reference vector REF is also
C     updated on exit.
C
C     first test for errors in the input arguments
C
      IERROR = 1
      NREC = 1
      IF (K.LT.1) THEN
         WRITE (P01REC,FMT=99999) K
      ELSE IF (LMAX.LT.2) THEN
         WRITE (P01REC,FMT=99998) LMAX
      ELSE IF (M.LE.0) THEN
         WRITE (P01REC,FMT=99997) M
      ELSE IF (M.GE.LMAX-MLAST) THEN
         NREC = 2
         WRITE (P01REC,FMT=99996) M, LMAX, MLAST
      ELSE IF (MLAST.LT.0) THEN
         WRITE (P01REC,FMT=99995) MLAST
      ELSE IF (IK.LT.K) THEN
         WRITE (P01REC,FMT=99994) IK, K
      ELSE
         IERROR = 0
         NREC = 0
C
C        test whether LREF is big enough
C
         JREF = (LMAX-1)*K*K + 2*K*LMAX + K
         IF (LREF.LT.JREF) THEN
            NREC = 2
            IERROR = 1
            WRITE (P01REC,FMT=99993) LREF, JREF
            GO TO 20
         END IF
C
C        call subroutine G13DKZ to update elements of PREDZ, SEFZ
C        and contents of REF
C
         LR1 = 1
         LR2 = LR1 + (LMAX-1)*K*K
         LR3 = LR2 + K*LMAX
         LR4 = LR3 + K*LMAX
C
         CALL G13DKZ(K,LMAX,Z,PREDZ,SEFZ,IK,M,REF(LR1),REF(LR2),REF(LR3)
     *               ,REF(LR4),WORK,MLAST,V,IERROR)
         MLAST = MLAST + M
C
         IF (IERROR.EQ.2) THEN
            NREC = 1
            WRITE (P01REC,FMT=99992)
         ELSE IF (IERROR.EQ.3) THEN
            NREC = 1
            WRITE (P01REC,FMT=99991)
         ELSE IF (IERROR.EQ.4) THEN
            NREC = 1
            WRITE (P01REC,FMT=99990)
         END IF
      END IF
   20 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT ('  ** On entry, K.lt.1 : K = ',I16)
99998 FORMAT ('  ** On entry, LMAX.lt.2 : LMAX = ',I16)
99997 FORMAT ('  ** On entry, M.le.0 : M = ',I16)
99996 FORMAT ('  ** On entry, M.ge.LMAX-MLAST : M = ',I16,'  and',/'  ',
     *       '   LMAX = ',I16,'  and MLAST = ',I16)
99995 FORMAT ('  ** On entry, MLAST.lt.0 : MLAST = ',I16)
99994 FORMAT ('  ** On entry, IK.lt.K : IK = ',I16,'  and K = ',I16)
99993 FORMAT ('  ** On entry, LREF is too small : LREF = ',I16,'  but ',
     *       'must be at',/'     least ',I16)
99992 FORMAT ('  ** On entry, some of the elements of the array REF ha',
     *       've been corrupted')
99991 FORMAT ('  ** On entry, one (or more) of the transformations req',
     *       'uested is invalid')
99990 FORMAT ('  ** The updated forecasts will overflow if computed')
      END
