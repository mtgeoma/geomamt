      SUBROUTINE E04NAU(MODE,UNITQ,QPHESS,N,NACTIV,NCTOTL,NFREE,NHESS,
     *                  NQ,NROWH,NCOLH,JADD,KACTIV,KFREE,ALFA,OBJQP,
     *                  GFIXED,GTP,CVEC,HESS,P,QTG,SCALE,X,ZY,WRK1,WRK2)
C     MARK 12 RE-ISSUE. NAG COPYRIGHT 1986.
C
C *********************************************************************
C     E04NAU  COMPUTES OR UPDATES...
C     (1)  OBJQP, THE VALUE OF THE QUADRATIC OBJECTIVE FUNCTION, AND
C (2)  THE VECTORS  Q(FREE)(T)G(FREE)  AND  G(FIXED),  WHERE  Q(FREE)
C       IS THE ORTHOGONAL FACTOR OF THE  A(FREE)  AND  A  IS THE MATRIX
C       OF CONSTRAINTS IN THE WORKING SET.  THESE VECTORS ARE STORED IN
C       ELEMENTS  1,2,...,NFREE  AND  NFREE+1,...,N,  RESPECTIVELY,  OF
C       THE ARRAY  QTG.
C     (3)  THE COMPONENT OF THE GRADIENT VECTOR CORRESPONDING TO A BOUND
C       CONSTRAINT THAT HAS JUST BEEN ADDED TO THE WORKING SET.
C
C     SYSTEMS OPTIMIZATION LABORATORY, STANFORD UNIVERSITY.
C     ORIGINAL VERSION OF OCTOBER 1982.
C *********************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, GFIXED, GTP, OBJQP
      INTEGER           JADD, MODE, N, NACTIV, NCOLH, NCTOTL, NFREE,
     *                  NHESS, NQ, NROWH
      LOGICAL           UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  CVEC(N), HESS(NROWH,NCOLH), P(N), QTG(N),
     *                  SCALE(NCTOTL), WRK1(N), WRK2(N), X(N), ZY(NQ,NQ)
      INTEGER           KACTIV(N), KFREE(N)
C     .. Subroutine Arguments ..
      EXTERNAL          QPHESS
C     .. Scalars in Common ..
      LOGICAL           SCLDQP
C     .. Local Scalars ..
      DOUBLE PRECISION  DELTAF, HALF, ONE, ZERO
      INTEGER           JTHCOL, NCOLZ
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          E04VDN, F06FBF, F06FCF, DAXPY, DCOPY
C     .. Common blocks ..
      COMMON            /JE04VC/SCLDQP
C     .. Data statements ..
      DATA              ZERO, HALF, ONE/0.0D+0, 0.5D+0, 1.0D+0/
C     .. Executable Statements ..
C
      JTHCOL = 0
      GO TO (20,40,60) MODE
C
C ---------------------------------------------------------------------
C     MODE = 1  ---  COMPUTE THE OBJECTIVE FUNCTION AND GRADIENT FROM
C                 SCRATCH.  ALLOW FOR A DIAGONAL SCALING OF  X.
C ---------------------------------------------------------------------
   20 CALL DCOPY(N,X,1,WRK1,1)
C
      IF (SCLDQP) CALL F06FCF(N,SCALE,1,WRK1,1)
C
      CALL QPHESS(N,NROWH,NCOLH,JTHCOL,HESS,WRK1,QTG)
      OBJQP = HALF*DDOT(N,QTG,1,WRK1,1) + DDOT(N,CVEC,1,WRK1,1)
      CALL DAXPY(N,ONE,CVEC,1,QTG,1)
C
      IF (SCLDQP) CALL F06FCF(N,SCALE,1,QTG,1)
C
C COMPUTE  Q(FREE)(T)(G(FREE)  AND  G(FIXED).  THE ELEMENTS OF  G(FREE)
C     ARE NOT STORED.
C
      CALL E04VDN(6,N,NACTIV,NCOLZ,NFREE,NQ,UNITQ,KACTIV,KFREE,QTG,ZY,
     *            WRK1)
C
      GO TO 80
C
C ---------------------------------------------------------------------
C MODE = 2  ---  IF THE QP OBJECTIVE FUNCTION IS REDUCED BY A POSITIVE
C                 STEP  ALFA,  OR  ALFA  IS NEGATIVE, UPDATE  OBJF,
C                 Q(FREE)(T)G(FREE)  AND  G(FIXED)  CORRESPONDING TO
C                 THE CHANGE,  X = X + ALFA P.
C ---------------------------------------------------------------------
   40 CALL QPHESS(N,NROWH,NCOLH,JTHCOL,HESS,P,WRK1)
C
      IF (SCLDQP) CALL F06FCF(N,SCALE,1,WRK1,1)
C
C     UPDATE  OBJQP.
C
      DELTAF = ALFA*GTP + HALF*ALFA*ALFA*DDOT(N,P,1,WRK1,1)
      IF (DELTAF.GT.ZERO .AND. ALFA.GT.ZERO) GO TO 100
      OBJQP = OBJQP + DELTAF
C
C     UPDATE THE ARRAY  QTG.  USE THE ARRAY  P  AS TEMPORARY WORK SPACE.
C
      CALL E04VDN(6,N,NACTIV,NCOLZ,NFREE,NQ,UNITQ,KACTIV,KFREE,WRK1,ZY,
     *            WRK2)
C
      CALL DAXPY(N,ALFA,WRK1,1,QTG,1)
      GO TO 80
C
C ---------------------------------------------------------------------
C MODE = 3  ---  COMPUTE THE  JADD-TH COMPONENT OF THE GRADIENT VECTOR.
C ---------------------------------------------------------------------
   60 JTHCOL = JADD
      CALL F06FBF(N,ZERO,WRK2,1)
      WRK2(JTHCOL) = ONE
      CALL QPHESS(N,NROWH,NCOLH,JTHCOL,HESS,WRK2,WRK1)
C
      IF (SCLDQP) CALL F06FCF(N,SCALE,1,WRK1,1)
      IF (SCLDQP) GFIXED = SCALE(JADD)*(DDOT(N,WRK1,1,X,1)+CVEC(JADD))
      IF ( .NOT. SCLDQP) GFIXED = DDOT(N,WRK1,1,X,1) + CVEC(JADD)
C
   80 NHESS = NHESS + 1
      RETURN
C
C     THE STEP  ALFA  DOES NOT DECREASE THE OBJECTIVE FUNCTION.
C
  100 MODE = -1
      RETURN
C
C     END OF E04NAU (QPGRAD)
      END
