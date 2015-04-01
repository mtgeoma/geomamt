      SUBROUTINE E04NCY(UNITQ,VERTEX,INFORM,K1,K2,NACTIV,NARTIF,NZ,
     *                  NFREE,NRANK,NREJTD,NRES,NGQ,N,LDZY,LDA,LDR,LDT,
     *                  ISTATE,KACTIV,KX,CONDMX,A,R,T,RES,GQ,ZY,W,C,S,
     *                  MSGLVL)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1075 (JUL 1993).
C
C     ******************************************************************
C     E04NCY  includes general constraints K1 thru K2 as new rows of
C     the TQ factorization stored in T, ZY.  If NRANK is nonzero, the
C     changes in Q are reflected in NRANK by N triangular factor R such
C     that
C                         C  =  P ( R ) Q,
C                                 ( 0 )
C     where  P  is orthogonal.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written  October-31-1984.
C     Level-2 matrix routines added 18-May-1988.
C     This version dated 27-April-1993.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CONDMX
      INTEGER           INFORM, K1, K2, LDA, LDR, LDT, LDZY, MSGLVL, N,
     *                  NACTIV, NARTIF, NFREE, NGQ, NRANK, NREJTD, NRES,
     *                  NZ
      LOGICAL           UNITQ, VERTEX
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(N), GQ(N,*), R(LDR,*), RES(N,*),
     *                  S(N), T(LDT,*), W(N), ZY(LDZY,*)
      INTEGER           ISTATE(*), KACTIV(N), KX(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, DTMAX, DTMIN
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION  CNDMAX, RNORM, ROWMAX, RTMAX
      INTEGER           I, IADD, IARTIF, IFIX, ISWAP, JADD, K, L, NZADD
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      EXTERNAL          DNRM2
C     .. External Subroutines ..
      EXTERNAL          E04NCV, F06FLF
C     .. Common blocks ..
      COMMON            /AX02ZA/WMACH
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
C     .. Save statement ..
      SAVE              /AX02ZA/
C     .. Executable Statements ..
      RTMAX = WMACH(8)
C
C     Estimate the condition number of the constraints that are not
C     to be refactorized.
C
      IF (NACTIV.EQ.0) THEN
         DTMAX = ZERO
         DTMIN = ONE
      ELSE
         CALL F06FLF(NACTIV,T(NACTIV,NZ+1),LDT-1,DTMAX,DTMIN)
      END IF
C
      DO 20 K = K1, K2
         IADD = KACTIV(K)
         JADD = N + IADD
         IF (NACTIV.LT.NFREE) THEN
C
            CALL E04NCV(UNITQ,INFORM,IFIX,IADD,JADD,NACTIV,NZ,NFREE,
     *                  NRANK,NRES,NGQ,N,LDA,LDZY,LDR,LDT,KX,CONDMX,A,R,
     *                  T,RES,GQ,ZY,W,C,S,MSGLVL)
C
            IF (INFORM.EQ.0) THEN
               NACTIV = NACTIV + 1
               NZ = NZ - 1
            ELSE
               ISTATE(JADD) = 0
               KACTIV(K) = -KACTIV(K)
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
C        If a vertex is required, add some temporary bounds.
C        We must accept the resulting condition number of the working
C        set.
C
         IF (VERTEX) THEN
            CNDMAX = RTMAX
            NZADD = NZ
            DO 80 IARTIF = 1, NZADD
               IF (UNITQ) THEN
                  IFIX = NFREE
                  JADD = KX(IFIX)
               ELSE
                  ROWMAX = ZERO
                  DO 60 I = 1, NFREE
                     RNORM = DNRM2(NZ,ZY(I,1),LDZY)
                     IF (ROWMAX.LT.RNORM) THEN
                        ROWMAX = RNORM
                        IFIX = I
                     END IF
   60             CONTINUE
                  JADD = KX(IFIX)
C
                  CALL E04NCV(UNITQ,INFORM,IFIX,IADD,JADD,NACTIV,NZ,
     *                        NFREE,NRANK,NRES,NGQ,N,LDA,LDZY,LDR,LDT,
     *                        KX,CNDMAX,A,R,T,RES,GQ,ZY,W,C,S,MSGLVL)
C
               END IF
               NFREE = NFREE - 1
               NZ = NZ - 1
               NARTIF = NARTIF + 1
               ISTATE(JADD) = 4
   80       CONTINUE
         END IF
      END IF
C
      NREJTD = K2 - NACTIV
C
      RETURN
C
C     End of  E04NCY. (LSADDS)
C
      END
