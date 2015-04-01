      SUBROUTINE H02BVF(N,M,A,LDA,BL,BU,X,CLAMDA,ISTATE,CRNAME,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='H02BVF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), BL(N+M), BU(N+M), CLAMDA(N+M), X(N)
      INTEGER           ISTATE(N+M)
      CHARACTER*8       CRNAME(N+M)
C     .. Scalars in Common ..
      DOUBLE PRECISION  BIGBND, BIGDX, BNDLOW, BNDUPP, TOLACT, TOLFEA,
     *                  TOLRNK
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLC(23), RPSVLC(30)
C     .. Local Scalars ..
      DOUBLE PRECISION  B1, B2, RES, RES2, V, WLAM
      INTEGER           IERR, IP, IS, J, K, NOUT, NPLIN, NREC
      CHARACTER*2       LS
C     .. Local Arrays ..
      CHARACTER*2       LSTATE(7)
      CHARACTER*80      REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      INTEGER           P01ABF
      EXTERNAL          DDOT, P01ABF
C     .. External Subroutines ..
      EXTERNAL          X04ABF, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Common blocks ..
      COMMON            /GE04MF/RPSVLC, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  TOLACT, TOLFEA, TOLRNK, RPADLC
C     .. Save statement ..
      SAVE              /GE04MF/
C     .. Data statements ..
      DATA              LSTATE(1)/'--'/, LSTATE(2)/'++'/
      DATA              LSTATE(3)/'FR'/, LSTATE(4)/'LL'/
      DATA              LSTATE(5)/'UL'/, LSTATE(6)/'EQ'/
      DATA              LSTATE(7)/'TF'/
C     .. Executable Statements ..
C
      CALL X04ABF(0,NOUT)
      IERR = 0
      NREC = 0
      IF (N.LT.1) THEN
         IERR = 1
         WRITE (REC,FMT=99999) N
         NREC = 2
         GO TO 40
      END IF
      IF (M.LE.0) THEN
         IERR = 1
         WRITE (REC,FMT=99998) M
         NREC = 2
         GO TO 40
      END IF
      IF (LDA.LT.MAX(1,M)) THEN
         IERR = 1
         WRITE (REC,FMT=99997) LDA, M
         NREC = 2
         GO TO 40
      END IF
C
      IF (BIGBND.EQ.-11111.0D0) BIGBND = 1.0D+20
C
      NPLIN = N + M
      WRITE (REC,FMT=99996)
      CALL X04BAY(NOUT,4,REC)
C
      DO 20 J = 1, NPLIN
         B1 = BL(J)
         B2 = BU(J)
         WLAM = CLAMDA(J)
         IS = ISTATE(J)
         LS = LSTATE(IS+3)
         IF (J.LE.N) THEN
C
C           Section 1 -- the variables  x.
C           ------------------------------
            K = J
            V = X(J)
C
         ELSE IF (J.LE.NPLIN) THEN
C
C           Section 2 -- the linear constraints  A*x.
C           -----------------------------------------
            IF (J.EQ.N+1) THEN
               WRITE (REC,FMT=99995)
               CALL X04BAY(NOUT,4,REC)
            END IF
C
            K = J - N
            V = DDOT(N,A(K,1),LDA,X,1)
         END IF
C
C        Print a line for the j-th variable or constraint.
C        -------------------------------------------------
         RES = V - B1
         RES2 = B2 - V
         IF (ABS(RES).GT.ABS(RES2)) RES = RES2
         IP = 1
         IF (B1.LE.(-BIGBND)) IP = 2
         IF (B2.GE.BIGBND) IP = IP + 2
C
         IF (IP.EQ.1) THEN
            WRITE (REC,FMT=99994) CRNAME(J), LS, V, B1, B2, WLAM, RES
         ELSE IF (IP.EQ.2) THEN
            WRITE (REC,FMT=99993) CRNAME(J), LS, V, B2, WLAM, RES
         ELSE IF (IP.EQ.3) THEN
            WRITE (REC,FMT=99992) CRNAME(J), LS, V, B1, WLAM, RES
         ELSE
            WRITE (REC,FMT=99991) CRNAME(J), LS, V, WLAM, RES
         END IF
         CALL X04BAF(NOUT,REC(1))
C
   20 CONTINUE
   40 IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
C
C     End of  H02BVF.
C
99999 FORMAT (' ** On entry, N.le.0:',/'    N = ',I16)
99998 FORMAT (' ** On entry, M.lt.0:',/'    M = ',I16)
99997 FORMAT (' ** On entry, LDA.lt.max(1,M):',/'    LDA = ',I16,'   M',
     *       ' = ',I16)
99996 FORMAT (//1X,'Varbl',3X,'State',5X,'Value',5X,'Lower Bound',3X,
     *       'Upper Bound',4X,'Lagr Mult',3X,'Residual',/)
99995 FORMAT (//1X,'L Con',3X,'State',5X,'Value',5X,'Lower Bound',3X,
     *       'Upper Bound',4X,'Lagr Mult',3X,'Residual',/)
99994 FORMAT (1X,A8,1X,A2,1X,1P,3G14.6,1P,2G12.4)
99993 FORMAT (1X,A8,1X,A2,1X,1P,G14.6,5X,'None',5X,1P,G14.6,1P,2G12.4)
99992 FORMAT (1X,A8,1X,A2,1X,1P,2G14.6,5X,'None',5X,1P,2G12.4)
99991 FORMAT (1X,A8,1X,A2,1X,1P,G14.6,5X,'None',10X,'None',5X,1P,2G12.4)
      END
