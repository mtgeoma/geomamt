      SUBROUTINE E04VDT(NULLR,UNITPG,UNITQ,N,NCLIN,NCLIN0,NCTOTL,NQ,
     *                  NROWA,NROWRT,NCOLRT,NCOLR,NCOLZ,NFREE,ISTATE,
     *                  KFREE,DINKY,GTP,PNORM,RDLAST,ZTGNRM,A,AP,P,QTG,
     *                  RT,V,ZY,WORK)
C     MARK 12 RE-ISSUE. NAG COPYRIGHT 1986.
C
C *********************************************************************
C     E04VDT COMPUTES THE FOLLOWING QUANTITIES FOR  E04MBY,  E04NAX  AND
C     LCCORE ...
C
C     (1) THE SEARCH DIRECTION  P  (AND ITS 2-NORM).
C     (2) THE VECTOR  V  SUCH THAT  R(T)V = - Z(T)G(FREE).  THIS VECTOR
C      IS REQUIRED BY  LCCORE  ONLY.
C     (3) THE VECTOR  AP,  WHERE  A  IS THE MATRIX OF LINEAR
C      CONSTRAINTS. AND, IF  NULLR  IS FALSE,
C     (4) THE  (NCOLR)-TH DIAGONAL ELEMENT OF THE CHOLESKY FACTOR OF THE
C      PROJECTED HESSIAN.
C
C     SYSTEMS OPTIMIZATION LABORATORY, STANFORD UNIVERSITY.
C     ORIGINAL VERSION OF DECEMBER 1982. REV. MAY 1983.
C *********************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DINKY, GTP, PNORM, RDLAST, ZTGNRM
      INTEGER           N, NCLIN, NCLIN0, NCOLR, NCOLRT, NCOLZ, NCTOTL,
     *                  NFREE, NQ, NROWA, NROWRT
      LOGICAL           NULLR, UNITPG, UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NROWA,N), AP(NCLIN0), P(N), QTG(N),
     *                  RT(NROWRT,NCOLRT), V(N), WORK(N), ZY(NQ,NQ)
      INTEGER           ISTATE(NCTOTL), KFREE(N)
C     .. Scalars in Common ..
      INTEGER           ISTART, MSG, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  ONE, ZERO
      INTEGER           I, IDIAG, J
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2
      EXTERNAL          DDOT, DNRM2
C     .. External Subroutines ..
      EXTERNAL          E04VDN, F04YAY, F06FBF, DAXPY, DCOPY, DSCAL,
     *                  X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Common blocks ..
      COMMON            /AE04VC/NOUT, MSG, ISTART
C     .. Data statements ..
      DATA              ZERO/0.0D+0/, ONE/1.0D+0/
C     .. Executable Statements ..
C
      CALL DCOPY(NCOLR,QTG,1,P,1)
      CALL DSCAL(NCOLR,(-ONE),P,1)
      IF (NULLR) GO TO 60
      RDLAST = RT(NCOLR,NCOLR)
C     ***
C     CORRECTION INSERTED BY MHW, 22 OCT 1985.
C     THIS ENSURES A NON-ZERO SEARCH DIRECTION.
C     ***
      IF (NCOLR.LT.NCOLZ .AND. ZTGNRM.LE.DINKY) P(NCOLR) = RDLAST
C
C ---------------------------------------------------------------------
C     SOLVE THE SYSTEM   R(T)R (PZ) = - Z(T)G(FREE).
C ---------------------------------------------------------------------
      IF (UNITPG) GO TO 20
C
C     PERFORM THE FORWARD SUBSTITUTION  R(T)V = - Z(T)G(FREE).
C
      IDIAG = 1
      CALL F04YAY(-1,NCOLR,RT,NROWRT,P,IDIAG)
      GO TO 40
C
C     THE PROJECTED GRADIENT IS A MULTIPLE OF THE UNIT VECTOR, THE
C     FORWARD SUBSTITUTION MAY BE AVOIDED.
C
   20 IF (ZTGNRM.LE.DINKY) P(NCOLR) = -ONE
      IF (ZTGNRM.GT.DINKY) P(NCOLR) = P(NCOLR)/RDLAST
C
C     PERFORM THE BACKWARD SUBSTITUTION   R(PZ) = P.
C
   40 CALL DCOPY(NCOLR,P,1,V,1)
      IDIAG = 1
      CALL F04YAY(1,NCOLR,RT,NROWRT,P,IDIAG)
C
C ---------------------------------------------------------------------
C     THE VECTOR  (PZ)  HAS BEEN COMPUTED.
C ---------------------------------------------------------------------
C     COMPUTE THE DIRECTIONAL DERIVATIVE  G(T)P = (GZ)(T)(PZ).
C
   60 GTP = DDOT(NCOLR,QTG,1,P,1)
C
C ---------------------------------------------------------------------
C     COMPUTE  P = Z * PZ.
C ---------------------------------------------------------------------
C     NACTIV  AND  KACTIV  ARE NOT USED IN  E04VDN.  N  AND  KFREE
C     SERVE AS ARGUMENTS FOR  NACTIV  AND  KACTIV.
C
      CALL E04VDN(1,N,N,NCOLR,NFREE,NQ,UNITQ,KFREE,KFREE,P,ZY,WORK)
C
      PNORM = DNRM2(NFREE,WORK,1)
      IF (MSG.GE.80) THEN
         WRITE (REC,FMT=99998)
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         DO 80 I = 1, N, 5
            WRITE (REC,FMT=99997) (P(J),J=I,MIN(N,I+4))
            CALL X04BAF(NOUT,REC(1))
   80    CONTINUE
      END IF
C
C ---------------------------------------------------------------------
C     COMPUTE  AP.
C ---------------------------------------------------------------------
      IF (NCLIN.EQ.0) GO TO 140
      CALL F06FBF(NCLIN,ZERO,AP,1)
      DO 100 J = 1, N
         IF (ISTATE(J).GT.0) GO TO 100
         CALL DAXPY(NCLIN,P(J),A(1,J),1,AP,1)
  100 CONTINUE
      IF (MSG.GE.80 .AND. NCLIN.GT.0) THEN
         WRITE (REC,FMT=99999)
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         DO 120 I = 1, NCLIN, 5
            WRITE (REC,FMT=99997) (AP(J),J=I,MIN(NCLIN,I+4))
            CALL X04BAF(NOUT,REC(1))
  120    CONTINUE
      END IF
C
  140 RETURN
C
C
C     END OF E04VDT (FINDP)
99999 FORMAT (/' //E04VDT//  AP ... ')
99998 FORMAT (/' //E04VDT//   P ... ')
99997 FORMAT (1P,5D15.5)
      END
