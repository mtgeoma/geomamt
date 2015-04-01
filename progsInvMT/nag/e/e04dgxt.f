      SUBROUTINE E04DGX(N,GAMMA,SJ,YJ,HJV,HJYJ,YJSJ,YJHYJ,VSJ,VHYJ,
     *                  HJP1V,ITER,DEBUG,IDBGCG)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 16A REVISED. IER-985 (JUN 1993).
C
C     E04DGX is a self scaled BFGS routine.
C
C     The formula here is that of page 16 of the CG paper,
C     except the k's there are J's here. Since no
C     approximate inverse Hessians are actually computed,
C     this routine returns the product of the updated
C     Hessian with a given vector.
C
C     For this reason, the calling routine must supply
C     the product of the given vector, say v, with
C     certain other vectors or matrices. For example
C     VSJ is the dot product of v with s(j), and VHYJ
C     is the product v(transpose)*h*y(j). Likewise,
C     HJV is the product of h and v ( h is old inverse
C     Hessian). On exit, HJP1V contains the product
C     of the updated Hessian inverse with v.
C
C     -- Written on 4-June-1986.
C     Sven Hammarling and Janet Welding, NAG Central Office.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  GAMMA, VHYJ, VSJ, YJHYJ, YJSJ
      INTEGER           IDBGCG, ITER, N
      LOGICAL           DEBUG
C     .. Array Arguments ..
      DOUBLE PRECISION  HJP1V(N), HJV(N), HJYJ(N), SJ(N), YJ(N)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  BETA, DELTA
      INTEGER           I, J
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
C     .. Executable Statements ..
C
      DELTA = (1.0D0+GAMMA*YJHYJ/YJSJ)*VSJ/YJSJ - GAMMA*VHYJ/YJSJ
      BETA = -GAMMA*VSJ/YJSJ
      DO 20 I = 1, N
         HJP1V(I) = GAMMA*HJV(I) + DELTA*SJ(I) + BETA*HJYJ(I)
   20 CONTINUE
      IF (DEBUG .AND. ITER.GE.IDBGCG) THEN
         WRITE (REC,FMT=99999)
         CALL X04BAY(NOUT,2,REC)
         DO 40 I = 1, N, 5
            WRITE (REC,FMT=99998) (HJP1V(J),J=I,MIN(I+4,N))
            CALL X04BAF(NOUT,REC(1))
   40    CONTINUE
      END IF
      RETURN
C
C     End of E04DGX (SSBFGS).
C
99999 FORMAT (/' //E04DGX// - N element vector HJP1V is')
99998 FORMAT (1X,5G15.7)
      END
