      SUBROUTINE E04NCT(UNITQ,N,NACTIV,NFREE,NRES,NGQ,NZ,NRZ,LDA,LDZY,
     *                  LDR,LDT,NRANK,JDEL,KDEL,KACTIV,KX,A,RES,R,T,GQ,
     *                  ZY,C,S)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 15B REVISED. IER-950 (NOV 1991).
C     MARK 16 REVISED. IER-1071 (JUL 1993).
C     MARK 17 REVISED. IER-1583 (JUN 1995).
C
C     ******************************************************************
C     E04NCT  updates the least-squares factor R and the factorization
C     A(free) (Z Y) = (0 T) when a regular, temporary or artificial
C     constraint is deleted from the working set.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October-1984.
C     Level-2 matrix routines added 25-Apr-1988.
C     This version of E04NCT dated 10-Sep-92.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           JDEL, KDEL, LDA, LDR, LDT, LDZY, N, NACTIV,
     *                  NFREE, NGQ, NRANK, NRES, NRZ, NZ
      LOGICAL           UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(N), GQ(N,*), R(LDR,*), RES(N,*),
     *                  S(N), T(LDT,*), ZY(LDZY,*)
      INTEGER           KACTIV(N), KX(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, DTMAX, DTMIN
C     .. Local Scalars ..
      DOUBLE PRECISION  CS, SN
      INTEGER           I, IR, ITDEL, JART, K, KA, LD, NPIV, NRZ1, NSUP,
     *                  NT
C     .. External Functions ..
      INTEGER           IDAMAX
      EXTERNAL          IDAMAX
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DSWAP, E04NBU, F06BAF, F06FBF, F06FLF,
     *                  F06QTF, F06QXF, F06QZZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Common blocks ..
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
C     .. Executable Statements ..
C
      IF (JDEL.GT.0) THEN
C        ---------------------------------------------------------------
C        Regular constraint or temporary bound deleted.
C        ---------------------------------------------------------------
C
         IF (JDEL.LE.N) THEN
C
C           Case 1.  A simple bound has been deleted.
C           =======  Columns NFREE+1 and IR of R must be swapped.
C
            IR = NZ + KDEL
C
            ITDEL = 1
            NFREE = NFREE + 1
C
            IF (NFREE.LT.IR) THEN
               KX(IR) = KX(NFREE)
               KX(NFREE) = JDEL
               IF (NRANK.GT.0) CALL E04NBU(N,NRES,NRANK,LDR,NFREE,IR,R,
     *                                     RES,C,S)
               CALL DSWAP(NGQ,GQ(NFREE,1),N,GQ(IR,1),N)
            END IF
C
            IF ( .NOT. UNITQ) THEN
C
C              Copy the incoming column of  A(free)  into the end of T.
C
               DO 20 KA = 1, NACTIV
                  I = KACTIV(KA)
                  T(KA,NFREE) = A(I,JDEL)
   20          CONTINUE
C
C              Expand Q by adding a unit row and column.
C
               IF (NFREE.GT.1) THEN
                  CALL F06FBF(NFREE-1,ZERO,ZY(NFREE,1),LDZY)
                  CALL F06FBF(NFREE-1,ZERO,ZY(1,NFREE),1)
               END IF
               ZY(NFREE,NFREE) = ONE
            END IF
         ELSE
C
C           Case 2.  A general constraint has been deleted.
C           =======
C
            ITDEL = KDEL
            NACTIV = NACTIV - 1
C
C           Delete row  kdel  of T and move up the ones below it.
C           T becomes reverse lower Hessenberg.
C
            DO 40 I = KDEL, NACTIV
               KACTIV(I) = KACTIV(I+1)
               LD = NFREE - I
               CALL DCOPY(I+1,T(I+1,LD),LDT,T(I,LD),LDT)
   40       CONTINUE
         END IF
C
         NZ = NZ + 1
C
         IF (NACTIV.EQ.0) THEN
            DTMAX = ONE
            DTMIN = ONE
         ELSE
C           ------------------------------------------------------------
C           Restore the NACTIV x (NACTIV+1) reverse-Hessenberg matrix  T
C           to reverse-triangular form.  The last NACTIV super-diagonal
C           elements are removed using a backward sweep of plane
C           rotations.  The rotation for the singleton in the first
C           column is generated separately.
C           ------------------------------------------------------------
            NSUP = NACTIV - ITDEL + 1
C
            IF (NSUP.GT.0) THEN
               NPIV = NFREE - ITDEL + 1
               IF (NSUP.GT.1) THEN
                  CALL DCOPY(NSUP-1,T(NACTIV-1,NZ+1),LDT-1,S(NZ+1),1)
                  CALL F06QZZ('Remove',NACTIV,1,NSUP,C(NZ+1),S(NZ+1),
     *                        T(1,NZ+1),LDT)
               END IF
C
               CALL F06BAF(T(NACTIV,NZ+1),T(NACTIV,NZ),CS,SN)
               T(NACTIV,NZ) = ZERO
               S(NZ) = -SN
               C(NZ) = CS
C
               CALL F06QXF('Right','Variable','Backwards',NFREE,NFREE,
     *                     NZ,NPIV,C,S,ZY,LDZY)
               CALL F06QXF('Left ','Variable','Backwards',NPIV,NGQ,NZ,
     *                     NPIV,C,S,GQ,N)
C
               NT = MIN(NRANK,NPIV)
C
               IF (NT.LT.NPIV .AND. NT.GT.0) THEN
C
C                 R is upper trapezoidal, pretend R is (NT x n) and
C                 apply the rotations in columns  max(NT,NZ)  thru NPIV.
C
                  CALL F06QXF('Right','Variable','Backwards',NT,N,
     *                        MAX(NT,NZ),NPIV,C,S,R,LDR)
               END IF
C
C              Apply the column transformations to the triangular part
C              of R.  A subdiagonal element is generated that must be
C              eliminated by a row rotation before the next column
C              transformation can be applied.
C
               IF (NZ.LT.NT) THEN
                  CALL F06QTF('Right',NT,NZ,NT,C,S,R,LDR)
               END IF
C
C              Apply the row rotations to the remaining rows of R.
C
               IF (N.GT.NT) CALL F06QXF('Left','Variable','Backwards',
     *                                  NT,N-NT,NZ,NT,C,S,R(1,NT+1),LDR)
C
               IF (NRES.GT.0) CALL F06QXF('Left','Variable','Backwards',
     *                                    NT,NRES,NZ,NT,C,S,RES,N)
C
            END IF
            CALL F06FLF(NACTIV,T(NACTIV,NZ+1),LDT-1,DTMAX,DTMIN)
         END IF
      END IF
C
      NRZ1 = NRZ + 1
C
      IF (NZ.GT.NRZ) THEN
         IF (JDEL.GT.0) THEN
            JART = NRZ1 - 1 + IDAMAX(NZ-NRZ1+1,GQ(NRZ1,1),1)
         ELSE
            JART = -JDEL
         END IF
C
         IF (JART.GT.NRZ1) THEN
C
C           Swap columns NRZ1 and JART of R.
C
            IF (UNITQ) THEN
               K = KX(NRZ1)
               KX(NRZ1) = KX(JART)
               KX(JART) = K
            ELSE
               CALL DSWAP(NFREE,ZY(1,NRZ1),1,ZY(1,JART),1)
            END IF
C
            CALL DSWAP(NGQ,GQ(NRZ1,1),N,GQ(JART,1),N)
            IF (NRANK.GT.0) CALL E04NBU(N,NRES,NRANK,LDR,NRZ1,JART,R,
     *                                  RES,C,S)
         END IF
      END IF
C
      NRZ = NRZ1
C
      RETURN
C
C
C     End of  E04NCT (LSDEL).
C
      END
