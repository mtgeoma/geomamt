      SUBROUTINE E04NBX(MSGLVL,NFREE,NROWA,N,NCLIN,NCTOTL,BIGBND,NAMED,
     *                  NAMES,NACTIV,ISTATE,KACTIV,KX,A,BL,BU,C,CLAMDA,
     *                  RLAMDA,X)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-591 (MAR 1988).
C     MARK 16 REVISED. IER-1058 (JUL 1993).
C
C     ******************************************************************
C     E04NBX   creates the expanded Lagrange multiplier vector CLAMDA.
C     If MSGLVL .EQ 1 or MSGLVL .GE. 10,  E04NBX prints  x,  A*x,
C     c(x),  their bounds, the multipliers, and the residuals (distance
C     to the nearer bound).
C
C     E04NBX is called by E04NCZ, E04UCZ and E04UPZ just before exiting.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original Fortran 77 version written  October 1984.
C     This version of  E04NBX  dated  30-Mar-1993.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LCMDBG
      PARAMETER         (LCMDBG=5)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND
      INTEGER           MSGLVL, N, NACTIV, NCLIN, NCTOTL, NFREE, NROWA
      LOGICAL           NAMED
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NROWA,*), BL(NCTOTL), BU(NCTOTL), C(*),
     *                  CLAMDA(NCTOTL), RLAMDA(N), X(N)
      INTEGER           ISTATE(NCTOTL), KACTIV(N), KX(N)
      CHARACTER*8       NAMES(*)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           CMDBG
C     .. Arrays in Common ..
      INTEGER           ICMDBG(LCMDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  B1, B2, RES, RES2, V, WLAM
      INTEGER           IP, IS, J, K, NFIXED, NPLIN, NZ
      CHARACTER*2       LS
      CHARACTER*5       ID3
      CHARACTER*8       ID4
C     .. Local Arrays ..
      CHARACTER*2       LSTATE(7)
      CHARACTER*5       ID(3)
      CHARACTER*80      REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          F06FBF, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /FE04NB/ICMDBG, CMDBG
C     .. Data statements ..
      DATA              ID(1)/'Varbl'/
      DATA              ID(2)/'L Con'/
      DATA              ID(3)/'N Con'/
      DATA              LSTATE(1)/'--'/, LSTATE(2)/'++'/
      DATA              LSTATE(3)/'FR'/, LSTATE(4)/'LL'/
      DATA              LSTATE(5)/'UL'/, LSTATE(6)/'EQ'/
      DATA              LSTATE(7)/'TF'/
C     .. Executable Statements ..
C
C
      NPLIN = N + NCLIN
      NZ = NFREE - NACTIV
C
C     Expand multipliers for bounds, linear and nonlinear constraints
C     into the  CLAMDA  array.
C
      CALL F06FBF(NCTOTL,ZERO,CLAMDA,1)
      NFIXED = N - NFREE
      DO 20 K = 1, NACTIV + NFIXED
         IF (K.LE.NACTIV) J = KACTIV(K) + N
         IF (K.GT.NACTIV) J = KX(NZ+K)
         CLAMDA(J) = RLAMDA(K)
   20 CONTINUE
C
      IF (MSGLVL.LT.10 .AND. MSGLVL.NE.1) RETURN
C
      WRITE (REC,FMT=99999)
      CALL X04BAY(IPRINT,4,REC)
      ID3 = ID(1)
C
      DO 40 J = 1, NCTOTL
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
               WRITE (REC,FMT=99998)
               CALL X04BAY(IPRINT,4,REC)
               ID3 = ID(2)
            END IF
C
            K = J - N
            V = DDOT(N,A(K,1),NROWA,X,1)
         ELSE
C
C           Section 3 -- the nonlinear constraints  c(x).
C           ---------------------------------------------
C
            IF (J.EQ.NPLIN+1) THEN
               WRITE (REC,FMT=99997)
               CALL X04BAY(IPRINT,4,REC)
               ID3 = ID(3)
            END IF
C
            K = J - NPLIN
            V = C(K)
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
         IF (NAMED) THEN
C
            ID4 = NAMES(J)
            IF (IP.EQ.1) THEN
               WRITE (REC,FMT=99996) ID4, LS, V, B1, B2, WLAM, RES
            ELSE IF (IP.EQ.2) THEN
               WRITE (REC,FMT=99995) ID4, LS, V, B2, WLAM, RES
            ELSE IF (IP.EQ.3) THEN
               WRITE (REC,FMT=99994) ID4, LS, V, B1, WLAM, RES
            ELSE
               WRITE (REC,FMT=99993) ID4, LS, V, WLAM, RES
            END IF
            CALL X04BAF(IPRINT,REC(1))
C
         ELSE
C
            IF (IP.EQ.1) THEN
               WRITE (REC,FMT=99992) ID3, K, LS, V, B1, B2, WLAM, RES
            ELSE IF (IP.EQ.2) THEN
               WRITE (REC,FMT=99991) ID3, K, LS, V, B2, WLAM, RES
            ELSE IF (IP.EQ.3) THEN
               WRITE (REC,FMT=99990) ID3, K, LS, V, B1, WLAM, RES
            ELSE
               WRITE (REC,FMT=99989) ID3, K, LS, V, WLAM, RES
            END IF
            CALL X04BAF(IPRINT,REC(1))
         END IF
   40 CONTINUE
      RETURN
C
C
C     End of  E04NBX. (CMPRT)
C
99999 FORMAT (//1X,'Varbl',1X,'State',5X,'Value',5X,'Lower Bound',3X,
     *       'Upper Bound',4X,'Lagr Mult',3X,'Residual',/)
99998 FORMAT (//1X,'L Con',1X,'State',5X,'Value',5X,'Lower Bound',3X,
     *       'Upper Bound',4X,'Lagr Mult',3X,'Residual',/)
99997 FORMAT (//1X,'N Con',1X,'State',5X,'Value',5X,'Lower Bound',3X,
     *       'Upper Bound',4X,'Lagr Mult',3X,'Residual',/)
99996 FORMAT (1X,A4,4X,A2,1X,1P,3G14.6,1P,2G12.4)
99995 FORMAT (1X,A4,4X,A2,1X,1P,G14.6,5X,'None',5X,1P,G14.6,1P,2G12.4)
99994 FORMAT (1X,A4,4X,A2,1X,1P,2G14.6,5X,'None',5X,1P,2G12.4)
99993 FORMAT (1X,A4,4X,A2,1X,1P,G14.6,5X,'None',10X,'None',5X,1P,2G12.4)
99992 FORMAT (1X,A1,I3,4X,A2,1X,1P,3G14.6,1P,2G12.4)
99991 FORMAT (1X,A1,I3,4X,A2,1X,1P,G14.6,5X,'None',5X,1P,G14.6,1P,
     *       2G12.4)
99990 FORMAT (1X,A1,I3,4X,A2,1X,1P,2G14.6,5X,'None',5X,2G12.4)
99989 FORMAT (1X,A1,I3,4X,A2,1X,1P,G14.6,5X,'None',10X,'None',5X,1P,
     *       2G12.4)
      END
