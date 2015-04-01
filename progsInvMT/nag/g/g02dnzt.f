      SUBROUTINE G02DNZ(IP,IRANK,B,COV,P,IMP,F,EST,STAT,SESTAT,T,TOL,WK,
     *                  IND)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     GIVES ESTIMATE AND SE OF ESTIMABLE FUNCTION
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SESTAT, STAT, T, TOL
      INTEGER           IMP, IND, IP, IRANK
      LOGICAL           EST
C     .. Array Arguments ..
      DOUBLE PRECISION  B(IP), COV((IP*IP+IP)/2), F(IP), P(*), WK(IP)
C     .. Local Scalars ..
      DOUBLE PRECISION  AT, D, TOLA
      INTEGER           IP2, LDP
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2, X02AJF
      EXTERNAL          DDOT, DNRM2, X02AJF
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DSPMV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      IF (IMP.EQ.0) THEN
         IP2 = 2*IP + 1
         LDP = IP
      ELSE
         IP2 = 1
         LDP = IMP
      END IF
C
C     CHECK IF FUNCTION IS ESTIMABLE
C
      IF (IP.EQ.IRANK) THEN
         EST = .TRUE.
      ELSE
         IF (TOL.LE.0.0D0) THEN
            TOLA = SQRT(X02AJF())
         ELSE
            TOLA = TOL
         END IF
         CALL DGEMV('N',IP-IRANK,IP,1.0D0,P(IP2+IRANK),LDP,F,1,0.0D0,WK,
     *              1)
         D = DNRM2(IP-IRANK,WK,1)
         IF (D.LT.TOLA) THEN
            EST = .TRUE.
         ELSE
            EST = .FALSE.
            IND = 0
         END IF
      END IF
      IF (EST) THEN
C
C        CALCULATE STATISTIC AND SE
C
         STAT = DDOT(IP,B,1,F,1)
         CALL DSPMV('U',IP,1.0D0,COV,F,1,0.0D0,WK,1)
         SESTAT = DDOT(IP,WK,1,F,1)
         IF (SESTAT.LE.0.0D0) THEN
            IND = 1
         ELSE
            IND = 0
            SESTAT = SQRT(SESTAT)
            T = STAT/SESTAT
            AT = ABS(T)
         END IF
      ELSE
         T = 0.0D0
         SESTAT = 0.0D0
         STAT = 0.0D0
      END IF
      END
