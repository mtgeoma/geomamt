      SUBROUTINE E04NFP(UNITQ,IT,N,NACTIV,NFREE,NGQ,NZ,NRZ,LDA,LDQ,LDT,
     *                  JDEL,KDEL,KACTIV,KX,A,T,GQM,Q,C,S)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1589 (JUN 1995).
C
C     ******************************************************************
C     E04NFP   updates the matrices  Z, Y, T, R  and  D  associated with
C     factorizations
C
C              A(free) * Q(free)  = (  0 T )
C                        Q(free)  = (  Z Y )
C
C     when a regular, temporary or artificial constraint is deleted
C     from the working set.
C
C     The  NACTIV x NACTIV  upper-triangular matrix  T  is stored
C     with its (1,1) element in position  (IT,JT)  of the array  T.
C
C     Original version written by PEG,  31-October-1984.
C     This version of  E04NFP  dated 14-Sep-1992.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           IT, JDEL, KDEL, LDA, LDQ, LDT, N, NACTIV, NFREE,
     *                  NGQ, NRZ, NZ
      LOGICAL           UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(N), GQM(N,*), Q(LDQ,*), S(N),
     *                  T(LDT,*)
      INTEGER           KACTIV(N), KX(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, DTMAX, DTMIN
C     .. Local Scalars ..
      DOUBLE PRECISION  CS, SN
      INTEGER           I, IR, ITDEL, J, JART, JT, K, L, NPIV, NRZ1,
     *                  NSUP
C     .. External Functions ..
      INTEGER           IDAMAX
      EXTERNAL          IDAMAX
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DSWAP, F06BAF, F06FBF, F06FLF, F06QRF,
     *                  F06QXF
C     .. Common blocks ..
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
C     .. Executable Statements ..
C
      JT = NZ + 1
C
      IF (JDEL.GT.0) THEN
C
C        Regular constraint or temporary bound deleted.
C
         IF (JDEL.LE.N) THEN
C
C           Case 1.  A simple bound has been deleted.
C           =======  Columns  NFREE+1  and  IR  of GQM' must be swapped.
C
            IR = NZ + KDEL
            ITDEL = NACTIV + 1
            NFREE = NFREE + 1
            IF (NFREE.LT.IR) THEN
               KX(IR) = KX(NFREE)
               KX(NFREE) = JDEL
               CALL DSWAP(NGQ,GQM(NFREE,1),N,GQM(IR,1),N)
            END IF
C
            IF ( .NOT. UNITQ) THEN
C
C              Copy the incoming column of  A(free)  into the end of  T.
C
               DO 20 K = 1, NACTIV
                  I = KACTIV(K)
                  T(NACTIV-K+1,NFREE) = A(I,JDEL)
   20          CONTINUE
C
C              Expand  Q  by adding a unit row and column.
C
               IF (NFREE.GT.1) THEN
                  CALL F06FBF(NFREE-1,ZERO,Q(NFREE,1),LDQ)
                  CALL F06FBF(NFREE-1,ZERO,Q(1,NFREE),1)
               END IF
               Q(NFREE,NFREE) = ONE
            END IF
         ELSE
C
C           Case 2.  A general constraint has been deleted.
C           =======
C
C           Delete row  ITDEL  of  T  and move up the ones below it.
C           T  becomes lower Hessenberg.
C
            ITDEL = KDEL
            DO 60 K = ITDEL, NACTIV
               J = JT + K - 1
               DO 40 L = ITDEL, K - 1
                  I = IT + L - 1
                  T(I,J) = T(I+1,J)
   40          CONTINUE
   60       CONTINUE
C
            DO 80 I = NACTIV - ITDEL + 1, NACTIV - 1
               KACTIV(I) = KACTIV(I+1)
   80       CONTINUE
            NACTIV = NACTIV - 1
         END IF
C
         NZ = NZ + 1
C
         IF (NACTIV.EQ.0) THEN
            DTMAX = ONE
            DTMIN = ONE
         ELSE
C           ------------------------------------------------------------
C           Restore the NACTIV x (NACTIV+1) upper-Hessenberg matrix  T
C           to upper-triangular form.  The  NSUP  super-diagonal
C           elements are removed by a backward sweep of rotations.
C           The rotation for the  (1,1)-th  element of  T  is generated
C           separately.
C           ------------------------------------------------------------
            NSUP = ITDEL - 1
C
            IF (NSUP.GT.0) THEN
               NPIV = JT + ITDEL - 1
               IF (NSUP.GT.1) THEN
                  CALL DCOPY(NSUP-1,T(IT+1,JT+1),LDT+1,S(JT+1),1)
                  CALL F06QRF('Right',NACTIV,1,NSUP,C(JT+1),S(JT+1),
     *                        T(IT,JT+1),LDT)
               END IF
C
               CALL F06BAF(T(IT,JT+1),T(IT,JT),CS,SN)
               T(IT,JT) = ZERO
               S(JT) = -SN
               C(JT) = CS
               CALL F06QXF('Right','Variable','Backwards',NFREE,NFREE,
     *                     NZ,NPIV,C,S,Q,LDQ)
               CALL F06QXF('Left ','Variable','Backwards',NPIV,NGQ,NZ,
     *                     NPIV,C,S,GQM,N)
            END IF
C
            JT = JT + 1
            CALL F06FLF(NACTIV,T(IT,JT),LDT+1,DTMAX,DTMIN)
         END IF
      END IF
C
      NRZ1 = NRZ + 1
C
      IF (NZ.GT.NRZ) THEN
         IF (JDEL.GT.0) THEN
            JART = NRZ1 - 1 + IDAMAX(NZ-NRZ1+1,GQM(NRZ1,1),1)
         ELSE
            JART = -JDEL
         END IF
C
         IF (JART.GT.NRZ1) THEN
C
C           Swap columns  NRZ1  and  JART  of  Q  and  GQM.
C
            IF (UNITQ) THEN
               K = KX(NRZ1)
               KX(NRZ1) = KX(JART)
               KX(JART) = K
            ELSE
               CALL DSWAP(NFREE,Q(1,NRZ1),1,Q(1,JART),1)
            END IF
C
            CALL DSWAP(NGQ,GQM(NRZ1,1),N,GQM(JART,1),N)
         END IF
      END IF
C
      NRZ = NRZ1
C
      RETURN
C
C     End of  E04NFP.  (RZDEL)
C
      END
