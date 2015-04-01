      SUBROUTINE E04NCL(HITCON,HITLOW,LINOBJ,UNITGZ,NCLIN,NRANK,NRZ,N,
     *                  LDR,JADD,NUMINF,ALFA,CTP,CTX,XNORM,AP,AX,BL,BU,
     *                  GQ,HZ,P,RES,R,X,WORK)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1064 (JUL 1993).
C     MARK 17 REVISED. IER-1576 (JUN 1995).
C
C     ******************************************************************
C     E04NCL  changes X to X + ALFA*P and updates CTX, AX, RES and GQ
C     accordingly.
C
C     If a bound was added to the working set,  move X exactly on to it,
C     except when a negative step was taken (E04UCG may have had to move
C     to some other closer constraint.)
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 27-December-1985.
C     This version dated 14-Sep-92.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, CTP, CTX, XNORM
      INTEGER           JADD, LDR, N, NCLIN, NRANK, NRZ, NUMINF
      LOGICAL           HITCON, HITLOW, LINOBJ, UNITGZ
C     .. Array Arguments ..
      DOUBLE PRECISION  AP(*), AX(*), BL(*), BU(*), GQ(*), HZ(*), P(N),
     *                  R(LDR,*), RES(*), WORK(*), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  BND
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      EXTERNAL          DNRM2
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DTRMV
C     .. Executable Statements ..
      CALL DAXPY(N,ALFA,P,1,X,1)
      IF (LINOBJ) CTX = CTX + ALFA*CTP
C
      IF (HITCON .AND. JADD.LE.N) THEN
         BND = BU(JADD)
         IF (HITLOW) BND = BL(JADD)
         IF (ALFA.GE.ZERO) X(JADD) = BND
      END IF
      XNORM = DNRM2(N,X,1)
C
      IF (NCLIN.GT.0) CALL DAXPY(NCLIN,ALFA,AP,1,AX,1)
C
      IF (NRZ.LE.NRANK) THEN
         IF (UNITGZ) THEN
            RES(NRZ) = RES(NRZ) - ALFA*HZ(NRZ)
         ELSE
            CALL DAXPY(NRZ,(-ALFA),HZ,1,RES,1)
         END IF
C
         IF (NUMINF.EQ.0) THEN
C
C           Update the transformed gradient GQ so that
C           GQ = GQ + ALFA*R'( HZ ).
C                            ( 0  )
C
            IF (UNITGZ) THEN
               CALL DAXPY(N-NRZ+1,ALFA*HZ(NRZ),R(NRZ,NRZ),LDR,GQ(NRZ),1)
            ELSE
               CALL DCOPY(NRZ,HZ,1,WORK,1)
               CALL DTRMV('U','T','N',NRZ,R,LDR,WORK,1)
               IF (NRZ.LT.N) CALL DGEMV('T',NRZ,N-NRZ,ONE,R(1,NRZ+1),
     *                                  LDR,HZ,1,ZERO,WORK(NRZ+1),1)
C
               CALL DAXPY(N,ALFA,WORK,1,GQ,1)
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of  E04NCL. (LSMOVE)
C
      END
