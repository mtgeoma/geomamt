      SUBROUTINE E04NCX(UNITQ,INFORM,NZ,NFREE,NRANK,NRES,NGQ,N,LDZY,LDA,
     *                  LDR,LDT,ISTATE,KX,CONDMX,A,R,T,RES,GQ,ZY,W,C,S,
     *                  MSGLVL)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1074 (JUL 1993).
C
C     ******************************************************************
C     E04NCX updates the factor R as KX is reordered to reflect the
C     status of the bound constraints given by ISTATE.  KX is reordered
C     so that the fixed variables come last.  One of two alternative
C     are used to reorder KX. One method needs fewer accesses to KX, the
C     other gives a matrix Rz with more rows and columns.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written  30-December-1985.
C     Level-2 matrix routines added 18-May-1988.
C     This version dated 27-April-1993.
C     ******************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CONDMX
      INTEGER           INFORM, LDA, LDR, LDT, LDZY, MSGLVL, N, NFREE,
     *                  NGQ, NRANK, NRES, NZ
      LOGICAL           UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(N), GQ(N,*), R(LDR,*), RES(N,*),
     *                  S(N), T(LDT,*), W(N), ZY(LDZY,*)
      INTEGER           ISTATE(*), KX(N)
C     .. Local Scalars ..
      INTEGER           IADD, IFIX, J, J2, JADD, K, L, LSTART, NACTV,
     *                  NFIXED
C     .. External Subroutines ..
      EXTERNAL          E04NBU, E04NCV
C     .. Executable Statements ..
      NFIXED = N - NFREE
C
      IF (NRANK.LT.N .AND. NRANK.GT.0) THEN
C        ---------------------------------------------------------------
C        R is specified but singular.  Try and keep the dimension of Rz
C        as large as possible.
C        ---------------------------------------------------------------
         NACTV = 0
         NFREE = N
         NZ = N
C
         J = N
C        +       WHILE (J .GT. 0  .AND.  N-NFREE .LT. NFIXED) DO
   20    IF (J.GT.0 .AND. N-NFREE.LT.NFIXED) THEN
            IF (ISTATE(J).GT.0) THEN
               JADD = J
               DO 40 IFIX = NFREE, 1, -1
                  IF (KX(IFIX).EQ.JADD) GO TO 60
   40          CONTINUE
C
C              Add bound JADD.
C
   60          CALL E04NCV(UNITQ,INFORM,IFIX,IADD,JADD,NACTV,NZ,NFREE,
     *                     NRANK,NRES,NGQ,N,LDA,LDZY,LDR,LDT,KX,CONDMX,
     *                     A,R,T,RES,GQ,ZY,W,C,S,MSGLVL)
C
               NFREE = NFREE - 1
               NZ = NZ - 1
            END IF
            J = J - 1
            GO TO 20
C           +       END WHILE
         END IF
      ELSE
C        ---------------------------------------------------------------
C        R is of full rank,  or is not specified.
C        ---------------------------------------------------------------
         IF (NFIXED.GT.0) THEN
C
C           Order KX so that the free variables come first.
C
            LSTART = NFREE + 1
            DO 120 K = 1, NFREE
               J = KX(K)
               IF (ISTATE(J).GT.0) THEN
                  DO 80 L = LSTART, N
                     J2 = KX(L)
                     IF (ISTATE(J2).EQ.0) GO TO 100
   80             CONTINUE
C
  100             KX(K) = J2
                  KX(L) = J
                  LSTART = L + 1
C
                  IF (NRANK.GT.0) CALL E04NBU(N,NRES,NRANK,LDR,K,L,R,
     *                                 RES,C,S)
               END IF
  120       CONTINUE
C
         END IF
         NZ = NFREE
      END IF
C
      RETURN
C
C     End of  E04NCX. (LSBNDS)
C
      END
