      SUBROUTINE E04VDU(MODFYG,ORTHOG,UNITQ,JDEL,KDEL,NACTIV,NCOLZ,
     *                  NFREE,N,NQ,NROWA,NROWRT,NCOLRT,KACTIV,KFREE,A,
     *                  QTG,RT,ZY)
C     MARK 12 RE-ISSUE. NAG COPYRIGHT 1986.
C
C *********************************************************************
C     E04VDU UPDATES THE FACTORIZATION OF THE MATRIX OF
C     CONSTRAINTS IN THE WORKING SET,  A(FREE) * (Z Y) = (0 T).
C
C     IF THERE ARE NO GENERAL CONSTRAINTS IN THE WORKING SET AND THE
C     MATRIX  Q = (Z Y)  IS THE IDENTITY,  Q  WILL NOT BE
C     TOUCHED.
C
C     SYSTEMS OPTIMIZATION LABORATORY, STANFORD UNIVERSITY.
C     VERSION OF DECEMBER 1981.  REV. OCT. 1982.
C *********************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           JDEL, KDEL, N, NACTIV, NCOLRT, NCOLZ, NFREE, NQ,
     *                  NROWA, NROWRT
      LOGICAL           MODFYG, ORTHOG, UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NROWA,N), QTG(N), RT(NROWRT,NCOLRT), ZY(NQ,NQ)
      INTEGER           KACTIV(N), KFREE(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, DTMAX, DTMIN
      INTEGER           ISTART, MSG, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  CS, ONE, SN, STORE, ZERO
      INTEGER           I, IBEGIN, IFREED, INCT, ISTORE, K, KA, KB, L,
     *                  LDIAG, LENQ, LENRT, NACTP1, NACTPI, NACTV1,
     *                  NCOLZ1, NFIXD1, NFREE1, NFREEI
C     .. Local Arrays ..
      CHARACTER*80      REC(4)
C     .. External Subroutines ..
      EXTERNAL          E04NAQ, E04NAR, F06FBF, F06FLF, DCOPY, X04BAF
C     .. Common blocks ..
      COMMON            /AE04VC/NOUT, MSG, ISTART
      COMMON            /HE04VC/ASIZE, DTMAX, DTMIN
C     .. Data statements ..
      DATA              ZERO/0.0D+0/, ONE/1.0D+0/
C     .. Executable Statements ..
C
      LENQ = NQ*(NQ-1) + 1
      IF (JDEL.GT.N) GO TO 80
C
C     ------------------------------------------------------------------
C     A SIMPLE BOUND IS BEING DELETED FROM THE WORKING SET.
C     ------------------------------------------------------------------
      IFREED = KDEL - NACTIV
      IF (MSG.GE.80) THEN
         WRITE (REC,FMT=99999) NACTIV, NCOLZ, NFREE, IFREED, JDEL, UNITQ
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         CALL X04BAF(NOUT,REC(3))
         CALL X04BAF(NOUT,REC(4))
      END IF
      NACTV1 = NACTIV
      NFREE1 = NFREE + 1
      IBEGIN = 1
      KFREE(NFREE1) = JDEL
C
C     ADD THE GRADIENT CORRESPONDING TO THE NEWLY-FREED VARIABLE TO THE
C     END OF  Q(FREE)(T)G(FREE).  THIS IS DONE BY INTERCHANGING THE
C     APPROPRIATE ELEMENTS OF  QTG  AND  KACTIV.
C
      IF ( .NOT. MODFYG) GO TO 20
      IF (IFREED.EQ.1) GO TO 20
      NFREEI = NFREE + IFREED
      NACTP1 = NACTIV + 1
      NACTPI = NACTIV + IFREED
      STORE = QTG(NFREE1)
      QTG(NFREE1) = QTG(NFREEI)
      QTG(NFREEI) = STORE
      ISTORE = KACTIV(NACTP1)
      KACTIV(NACTP1) = KACTIV(NACTPI)
      KACTIV(NACTPI) = ISTORE
C
C     COPY THE INCOMING COLUMN OF  A  INTO THE END OF  T.
C
   20 IF (UNITQ) GO TO 120
      IF (NACTIV.EQ.0) GO TO 60
C
      DO 40 KA = 1, NACTIV
         I = KACTIV(KA)
         RT(KA,NFREE1) = A(I,JDEL)
   40 CONTINUE
C
C     EXPAND  Q  BY ADDING A UNIT ROW AND COLUMN.
C
   60 CALL F06FBF(NFREE,ZERO,ZY(NFREE1,1),NQ)
      CALL F06FBF(NFREE,ZERO,ZY(1,NFREE1),1)
      ZY(NFREE1,NFREE1) = ONE
      GO TO 120
C
C     ------------------------------------------------------------------
C     A GENERAL CONSTRAINT IS BEING DELETED FROM THE WORKING SET.
C     ------------------------------------------------------------------
   80 IF (MSG.GE.80) THEN
         WRITE (REC,FMT=99998) NACTIV, NCOLZ, NFREE, KDEL, JDEL, UNITQ
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         CALL X04BAF(NOUT,REC(3))
         CALL X04BAF(NOUT,REC(4))
      END IF
      NACTV1 = NACTIV - 1
      NFREE1 = NFREE
      IBEGIN = KDEL
      IF (KDEL.GT.NACTV1) GO TO 120
C
C     DELETE A ROW OF  T  AND MOVE THE ONES BELOW IT UP.
C
      DO 100 I = KDEL, NACTV1
         KACTIV(I) = KACTIV(I+1)
         LENRT = NROWRT*I + 1
         LDIAG = NFREE - I
         CALL DCOPY(I+1,RT(I+1,LDIAG),NROWRT,RT(I,LDIAG),NROWRT)
  100 CONTINUE
C
C     ------------------------------------------------------------------
C     ELIMINATE THE SUPER-DIAGONAL ELEMENTS OF  T,
C     USING A BACKWARD SWEEP OF 2*2 TRANFORMATIONS.
C     ------------------------------------------------------------------
  120 IF (IBEGIN.GT.NACTV1) GO TO 160
      K = NFREE1 - IBEGIN
      L = NACTV1 - IBEGIN
C
      DO 140 I = IBEGIN, NACTV1
         CALL E04NAQ(ORTHOG,RT(I,K+1),RT(I,K),CS,SN)
         IF (L.GT.0) CALL E04NAR(ORTHOG,L,RT(I+1,K+1),L,1,RT(I+1,K),L,1,
     *                           CS,SN)
         IF (NACTV1.GT.0) CALL E04NAR(ORTHOG,NFREE1,ZY(1,K+1),NQ,1,ZY(1,
     *                                K),NQ,1,CS,SN)
         IF (MODFYG) CALL E04NAR(ORTHOG,1,QTG(K+1),1,1,QTG(K),1,1,CS,SN)
         K = K - 1
         L = L - 1
  140 CONTINUE
C
C     ------------------------------------------------------------------
C     COMPRESS THE ELEMENTS OF  KACTIV  CORRESPONDING TO FIXED
C     VARIABLES.
C     ------------------------------------------------------------------
  160 NFIXD1 = N - NFREE1
      KB = NACTV1 + 1
      IF (NFIXD1.EQ.0) GO TO 200
      DO 180 K = 1, NFIXD1
         KACTIV(KB) = KACTIV(KB+1)
         KB = KB + 1
  180 CONTINUE
C
C     ------------------------------------------------------------------
C     ESTIMATE THE CONDITION NUMBER OF  T.
C     ------------------------------------------------------------------
  200 NCOLZ1 = NCOLZ + 1
      LENRT = NROWRT*(NACTV1-1) + 1
      INCT = NROWRT - 1
      IF (NACTV1.GT.0) CALL F06FLF(NACTV1,RT(NACTV1,NCOLZ1+1),INCT,
     *                             DTMAX,DTMIN)
C
      RETURN
C
C
C     END OF E04VDU (DELCON)
99999 FORMAT (/' //E04VDU//  SIMPLE BOUND DELETED.',/' //E04VDU//  NAC',
     *  'TIV NCOLZ NFREE IFREED JDEL UNITQ',/' //E04VDU//  ',3I6,I7,I5,
     *  L6)
99998 FORMAT (/' //E04VDU//  GENERAL CONSTRAINT DELETED.',/' //E04VDU/',
     *  '/  NACTIV NCOLZ NFREE  KDEL  JDEL UNITQ',/' //E04VDU//  ',5I6,
     *  L6)
      END
