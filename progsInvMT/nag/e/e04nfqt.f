      SUBROUTINE E04NFQ(UNITQ,VERTEX,K1,K2,IT,NACTIV,NARTIF,NZ,NFREE,
     *                  NREJTD,NGQ,N,LDQ,LDA,LDT,ISTATE,KACTIV,KX,
     *                  CONDMX,A,T,GQM,Q,W,C,S,MSGLVL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1590 (JUN 1995).
C
C     ******************************************************************
C     E04NFQ  includes general constraints  K1  thru  K2  as new rows of
C     the  TQ  factorization:
C              A(free) * Q(free)  = (  0 T )
C                        Q(free)  = (  Z Y )
C
C     a) The  NACTIV x NACTIV  upper-triangular matrix  T  is stored
C        with its (1,1) element in position  (IT,JT)  of the array  T.
C
C     Original version written by PEG,  October-31-1984.
C     This version of  E04NFQ  dated  7-Jul-1989.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CONDMX
      INTEGER           IT, K1, K2, LDA, LDQ, LDT, MSGLVL, N, NACTIV,
     *                  NARTIF, NFREE, NGQ, NREJTD, NZ
      LOGICAL           UNITQ, VERTEX
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(N), GQM(N,*), Q(LDQ,*), S(N),
     *                  T(LDT,*), W(N)
      INTEGER           ISTATE(*), KACTIV(N), KX(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, DTMAX, DTMIN
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION  CNDMAX, COND, DELTA, DRZZ, DTNEW, RNORM, ROWMAX,
     *                  RTMAX, TDTMAX, TDTMIN
      INTEGER           I, IADD, IARTIF, IFIX, INFORM, ISWAP, J, JADD,
     *                  JT, K, L, NZADD
      LOGICAL           OVERFL, RSET
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, F06BLF
      EXTERNAL          DNRM2, F06BLF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, DGER, E04NBW, E04NFR, F06FLF,
     *                  F06FRF, F06QHF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Common blocks ..
      COMMON            /AX02ZA/WMACH
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
C     .. Save statement ..
      SAVE              /AX02ZA/
C     .. Executable Statements ..
C
      RTMAX = WMACH(8)
C
      JT = NZ + 1
C
C     Estimate the condition number of the constraints already
C     factorized.
C
      IF (NACTIV.EQ.0) THEN
         DTMAX = ZERO
         DTMIN = ONE
         IF (UNITQ) THEN
C
C           First general constraint added.  Set  Q = I.
C
            CALL F06QHF('General',NFREE,NFREE,ZERO,ONE,Q,LDQ)
            UNITQ = .FALSE.
         END IF
      ELSE
         CALL F06FLF(NACTIV,T(IT,JT),LDT+1,DTMAX,DTMIN)
      END IF
C
      DO 20 K = K1, K2
         IADD = KACTIV(K)
         JADD = N + IADD
         IF (NACTIV.LT.NFREE) THEN
C
            OVERFL = .FALSE.
C
C           Transform the incoming row of  A  by  Q'.
C
            CALL DCOPY(N,A(IADD,1),LDA,W,1)
            CALL E04NBW(8,N,NZ,NFREE,LDQ,UNITQ,KX,W,Q,S)
C
C           Check that the incoming row is not dependent upon those
C           already in the working set.
C
            DTNEW = DNRM2(NZ,W,1)
            IF (NACTIV.EQ.0) THEN
C
C              This is the first general constraint in the working set.
C
               COND = F06BLF(ASIZE,DTNEW,OVERFL)
               TDTMAX = DTNEW
               TDTMIN = DTNEW
            ELSE
C
C              There are already some general constraints in the working
C              set. Update the estimate of the condition number.
C
               TDTMAX = MAX(DTNEW,DTMAX)
               TDTMIN = MIN(DTNEW,DTMIN)
               COND = F06BLF(TDTMAX,TDTMIN,OVERFL)
            END IF
C
            IF (COND.GE.CONDMX .OR. OVERFL) THEN
C              ---------------------------------------------------------
C              This constraint appears to be dependent on those already
C              in the working set.  Skip it.
C              ---------------------------------------------------------
               ISTATE(JADD) = 0
               KACTIV(K) = -KACTIV(K)
            ELSE
               IF (NZ.GT.1) THEN
C                 ------------------------------------------------------
C                 Use a single column transformation to reduce the first
C                 NZ-1  elements of  W  to zero.
C                 ------------------------------------------------------
C                 Apply the Householder reflection  I  -  W W'.
C                 The reflection is applied to  Z  and GQM so that
C                    Y  =    Z  * W,   Z    =  Z    -  Y W'  and
C                    Y  =  GQM' * W,   GQM  =  GQM  -  W Y',
C                 where  W = WRK1 (from Householder),
C                 and    Y = WRK2 (workspace).
C
C                 Note that DELTA  has to be stored after the reflection
C                 is used.
C
                  DELTA = W(NZ)
                  CALL F06FRF(NZ-1,DELTA,W,1,ZERO,W(NZ))
                  IF (W(NZ).GT.ZERO) THEN
C
                     CALL DGEMV('N',NFREE,NZ,ONE,Q,LDQ,W,1,ZERO,S,1)
                     CALL DGER(NFREE,NZ,(-ONE),S,1,W,1,Q,LDQ)
C
                     IF (NGQ.GT.0) THEN
                        CALL DGEMV('T',NZ,NGQ,ONE,GQM,N,W,1,ZERO,S,1)
                        CALL DGER(NZ,NGQ,(-ONE),W,1,S,1,GQM,N)
                     END IF
                  END IF
C
                  W(NZ) = DELTA
               END IF
               IT = IT - 1
               JT = JT - 1
               NACTIV = NACTIV + 1
               NZ = NZ - 1
               CALL DCOPY(NACTIV,W(JT),1,T(IT,JT),LDT)
               DTMAX = TDTMAX
               DTMIN = TDTMIN
            END IF
         END IF
   20 CONTINUE
C
      IF (NACTIV.LT.K2) THEN
C
C        Some of the constraints were classed as dependent and not
C        included in the factorization.  Re-order the part of  KACTIV
C        that holds the indices of the general constraints in the
C        working set.  Move accepted indices to the front and shift
C        rejected indices (with negative values) to the end.
C
         L = K1 - 1
         DO 40 K = K1, K2
            I = KACTIV(K)
            IF (I.GE.0) THEN
               L = L + 1
               IF (L.NE.K) THEN
                  ISWAP = KACTIV(L)
                  KACTIV(L) = I
                  KACTIV(K) = ISWAP
               END IF
            END IF
   40    CONTINUE
C
C        If a vertex is required,  add some temporary bounds.
C        We must accept the resulting condition number of the working
C        set.
C
         IF (VERTEX) THEN
            RSET = .FALSE.
            CNDMAX = RTMAX
            DRZZ = ONE
            NZADD = NZ
            DO 80 IARTIF = 1, NZADD
               IF (UNITQ) THEN
                  IFIX = NFREE
                  JADD = KX(IFIX)
               ELSE
                  ROWMAX = ZERO
                  DO 60 I = 1, NFREE
                     RNORM = DNRM2(NZ,Q(I,1),LDQ)
                     IF (ROWMAX.LT.RNORM) THEN
                        ROWMAX = RNORM
                        IFIX = I
                     END IF
   60             CONTINUE
                  JADD = KX(IFIX)
C
                  CALL E04NFR(UNITQ,RSET,INFORM,IFIX,IADD,JADD,IT,
     *                        NACTIV,NZ,NFREE,NZ,NGQ,N,LDA,LDQ,LDT,LDT,
     *                        KX,CNDMAX,DRZZ,A,T,T,GQM,Q,W,C,S,MSGLVL)
               END IF
               NFREE = NFREE - 1
               NZ = NZ - 1
               NARTIF = NARTIF + 1
               ISTATE(JADD) = 4
   80       CONTINUE
         END IF
C
         IF (IT.GT.1) THEN
C           ------------------------------------------------------------
C           If some dependent constraints were rejected,  move  T  to
C           the top of the array  T.
C           ------------------------------------------------------------
            DO 120 K = 1, NACTIV
               J = NZ + K
               DO 100 I = 1, K
                  T(I,J) = T(IT+I-1,J)
  100          CONTINUE
  120       CONTINUE
         END IF
      END IF
C
      NREJTD = K2 - NACTIV
C
      RETURN
C
C     End of  E04NFQ.  (RZADDS)
C
      END
