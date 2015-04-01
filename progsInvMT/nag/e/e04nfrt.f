      SUBROUTINE E04NFR(UNITQ,RSET,INFORM,IFIX,IADD,JADD,IT,NACTIV,NZ,
     *                  NFREE,NRZ,NGQ,N,LDA,LDQ,LDR,LDT,KX,CONDMX,DRZZ,
     *                  A,R,T,GQM,Q,W,C,S,MSGLVL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1591 (JUN 1995).
C
C     ******************************************************************
C     E04NFR  updates the matrices  Z, Y, T, R  and  D  associated with
C     factorizations
C
C              A(free) * Q(free)  = (  0 T )
C                        Q(free)  = (  Z Y )
C                      R' *D * R  =   Hz
C
C     a) The matrices  R  and  T  are upper triangular.
C     b) The arrays  T  and  R  may be the same array.
C     c) The  NACTIV x NACTIV  upper-triangular matrix  T  is stored
C        with its (1,1) element in position  (IT,JT) of the
C        array  T.   The integer  JT  is always  NZ+1.  During regular
C        changes to the working set,  IT = 1;  when several constraints
C        are added simultaneously,  IT  points to the first row of the
C        existing  T.
C     d) The matrix  R  is stored in the first  NZ x NZ  rows
C        and columns of the  NFREE x NFREE  leading principal minor of
C        the array  R.
C     e) If  RSET  is  false,   R  is not touched.
C
C     There are three separate cases to consider (although each case
C     shares code with another)...
C
C     (1) A free variable becomes fixed on one of its bounds when there
C         are already some general constraints in the working set.
C
C     (2) A free variable becomes fixed on one of its bounds when there
C         are only bound constraints in the working set.
C
C     (3) A general constraint (corresponding to row  IADD  of  A) is
C         added to the working set.
C
C     In cases (1) and (2), we assume that  KX(IFIX) = JADD.
C     In all cases,  JADD  is the index of the constraint being added.
C
C     If there are no general constraints in the working set,  the
C     matrix  Q = (Z Y)  is the identity and will not be touched.
C
C     If  NGQ .gt. 0,  the column transformations are applied to the
C     columns of the  (NGQ x N)  matrix  GQM'.
C
C     Original version written by PEG,  31-October-1984.
C     This version of  E04NFR  dated  7-Jul-1989.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CONDMX, DRZZ
      INTEGER           IADD, IFIX, INFORM, IT, JADD, LDA, LDQ, LDR,
     *                  LDT, MSGLVL, N, NACTIV, NFREE, NGQ, NRZ, NZ
      LOGICAL           RSET, UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(N), GQM(N,*), Q(LDQ,*), R(LDR,*),
     *                  S(N), T(LDT,*), W(N)
      INTEGER           KX(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, DTMAX, DTMIN, EPSPT3, EPSPT5, EPSPT8,
     *                  EPSPT9
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  COND, CONDBD, DTNEW, TDTMAX, TDTMIN
      INTEGER           I, J, JT, K, NANEW, NPIV, NSUP
      LOGICAL           BOUND, OVERFL
C     .. Local Arrays ..
      CHARACTER*80      REC(5)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, F06BLF
      EXTERNAL          DNRM2, F06BLF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DSCAL, E04NBW, E04NFM, F06FLF, F06FQF,
     *                  F06QHF, F06QKF, F06QRF, F06QVF, F06QXF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
C     .. Executable Statements ..
C
C     If the condition estimator of the updated T is greater than
C     CONDBD,  a warning message is printed.
C
      CONDBD = ONE/EPSPT9
C
      OVERFL = .FALSE.
      BOUND = JADD .LE. N
      JT = NZ + 1
C
      IF (BOUND) THEN
C        ===============================================================
C        A simple bound has entered the working set.  IADD is not used.
C        ===============================================================
         NANEW = NACTIV
C
         IF (UNITQ) THEN
C
C           Q is not stored, but  KX  defines an ordering of the columns
C           of the identity matrix that implicitly define Q.
C           Define the sequence of pairwise interchanges P that moves
C           the newly-fixed variable to position  NFREE.
C           Reorder  KX  accordingly.
C
            DO 20 I = 1, NFREE - 1
               IF (I.GE.IFIX) THEN
                  W(I) = I + 1
                  KX(I) = KX(I+1)
               ELSE
                  W(I) = I
               END IF
   20       CONTINUE
         ELSE
C           ------------------------------------------------------------
C           Q  is stored explicitly.
C           ------------------------------------------------------------
C           Set  W = the  (IFIX)-th  row of  Q.
C           Move the  (NFREE)-th  row of  Q  to position IFIX.
C
            CALL DCOPY(NFREE,Q(IFIX,1),LDQ,W,1)
            IF (IFIX.LT.NFREE) THEN
               CALL DCOPY(NFREE,Q(NFREE,1),LDQ,Q(IFIX,1),LDQ)
               KX(IFIX) = KX(NFREE)
            END IF
         END IF
         KX(NFREE) = JADD
      ELSE
C        ===============================================================
C        A general constraint has entered the working set.
C        IFIX is not used.
C        ===============================================================
         NANEW = NACTIV + 1
C
C        Transform the incoming row of A by Q'.
C
         CALL DCOPY(N,A(IADD,1),LDA,W,1)
         CALL E04NBW(8,N,NZ,NFREE,LDQ,UNITQ,KX,W,Q,C)
C
C        Check that the incoming row is not dependent upon those
C        already in the working set.
C
         DTNEW = DNRM2(NZ,W,1)
         IF (NACTIV.EQ.0) THEN
C
C           This is the only general constraint in the working set.
C
            COND = F06BLF(ASIZE,DTNEW,OVERFL)
            TDTMAX = DTNEW
            TDTMIN = DTNEW
         ELSE
C
C           There are already some general constraints in the working
C           set.  Update the estimate of the condition number.
C
            TDTMAX = MAX(DTNEW,DTMAX)
            TDTMIN = MIN(DTNEW,DTMIN)
            COND = F06BLF(TDTMAX,TDTMIN,OVERFL)
         END IF
C
         IF (COND.GT.CONDMX .OR. OVERFL) GO TO 80
C
         IF (UNITQ) THEN
C
C           First general constraint added.  Set  Q = I.
C
            CALL F06QHF('General',NFREE,NFREE,ZERO,ONE,Q,LDQ)
            UNITQ = .FALSE.
            IT = 0
         END IF
      END IF
C
      IF (BOUND) THEN
         NPIV = NFREE
      ELSE
         NPIV = NZ
      END IF
C
      IF (UNITQ) THEN
C        ---------------------------------------------------------------
C        The orthogonal matrix  Q  (i.e.,  Q) is not stored explicitly.
C        Apply  P, the sequence of pairwise interchanges that moves the
C        newly-fixed variable to position  NFREE.
C        ---------------------------------------------------------------
         IF (NGQ.GT.0) CALL F06QKF('Left','Transpose',NFREE-1,W,NGQ,GQM,
     *                             N)
C
         IF (RSET) THEN
C
C           Apply the pairwise interchanges to  Rz.
C           The subdiagonal elements generated by this process are
C           stored in  S(IFIX), S(2), ..., S(NRZ-1).
C
            NSUP = NRZ - IFIX
            CALL E04NFM('Right',NRZ,IFIX,NRZ,S,R,LDR)
         END IF
      ELSE
C        ---------------------------------------------------------------
C        The matrix  Q  is stored explicitly.
C        Define a sweep of plane rotations P such that
C                           PW = beta*e(NPIV).
C        The rotations are applied in the planes (1, 2), (2, 3), ...,
C        (NPIV-1, NPIV).  The rotations must be applied to Q, GQM', R
C        and T.
C        ---------------------------------------------------------------
         CALL F06FQF('Varble','Forwrds',NPIV-1,W(NPIV),W,1,C,S)
C
         IF (NGQ.GT.0) CALL F06QXF('Left ','Variable','Forwards',NPIV,
     *                             NGQ,1,NPIV,C,S,GQM,N)
         CALL F06QXF('Right','Variable','Forwards',NFREE,NFREE,1,NPIV,C,
     *               S,Q,LDQ)
C
         IF (RSET) THEN
C
C           Apply the rotations to the triangular part of R.
C           The subdiagonal elements generated by this process are
C           stored in  S(1),  S(2), ..., S(NRZ-1).
C
            NSUP = NRZ - 1
            CALL F06QVF('Right',NRZ,1,NRZ,C,S,R,LDR)
         END IF
      END IF
C
      IF (RSET) THEN
C        ---------------------------------------------------------------
C        Eliminate the  NSUP  subdiagonal elements of  R  stored in
C        S(NRZ-NSUP), ..., S(NRZ-1)  with a left-hand sweep of rotations
C        in planes (NRZ-NSUP, NRZ-NSUP+1), ..., (NRZ-1, NRZ).
C        ---------------------------------------------------------------
         CALL F06QRF('Left ',NRZ,NRZ-NSUP,NRZ,C,S,R,LDR)
C
         IF (NSUP.GT.0 .AND. DRZZ.NE.ONE) THEN
            DRZZ = C(NRZ-1)**2 + DRZZ*S(NRZ-1)**2
         END IF
      END IF
C
      IF ( .NOT. UNITQ) THEN
         IF (BOUND) THEN
C           ------------------------------------------------------------
C           Bound constraint added.   The rotations affect columns
C           NZ+1  thru  NFREE  of  GQM'  and  T.
C           ------------------------------------------------------------
C           The last row and column of  Q  has been transformed to plus
C           or minus the unit vector  e(NFREE).  We can reconstitute the
C           column of GQM' corresponding to the new fixed variable.
C
            IF (W(NFREE).LT.ZERO) THEN
               IF (NGQ.GT.0) CALL DSCAL(NGQ,-ONE,GQM(NFREE,1),N)
            END IF
C
            IF (NACTIV.GT.0) THEN
               T(IT,JT-1) = S(JT-1)*T(IT,JT)
               T(IT,JT) = C(JT-1)*T(IT,JT)
C
               IF (NACTIV.GT.1) THEN
                  CALL F06QVF('Right',NACTIV,1,NACTIV,C(JT),S(JT),
     *                        T(IT,JT),LDT)
                  CALL DCOPY(NACTIV-1,S(JT),1,T(IT+1,JT),LDT+1)
               END IF
C
               JT = JT - 1
               CALL F06FLF(NACTIV,T(IT,JT),LDT+1,TDTMAX,TDTMIN)
               COND = F06BLF(TDTMAX,TDTMIN,OVERFL)
            END IF
         ELSE
C           ------------------------------------------------------------
C           General constraint added.  Install  W  at the front of  T.
C           If there is no room,  shift all the rows down one position.
C           ------------------------------------------------------------
            IT = IT - 1
            IF (IT.LE.0) THEN
               IT = 1
               DO 60 K = 1, NACTIV
                  J = JT + K - 1
                  DO 40 I = K, 1, -1
                     T(I+1,J) = T(I,J)
   40             CONTINUE
   60          CONTINUE
            END IF
            JT = JT - 1
            CALL DCOPY(NANEW,W(JT),1,T(IT,JT),LDT)
         END IF
      END IF
C
C     ==================================================================
C     Prepare to exit.  Check the magnitude of the condition estimator.
C     ==================================================================
   80 IF (NANEW.GT.0) THEN
         IF (COND.LT.CONDMX .AND. .NOT. OVERFL) THEN
C
C           The factorization has been successfully updated.
C
            INFORM = 0
            DTMAX = TDTMAX
            DTMIN = TDTMIN
            IF (COND.GE.CONDBD) THEN
               IF (MSGLVL.GT.0) THEN
                  WRITE (REC,FMT=99999) JADD
                  CALL X04BAY(IPRINT,5,REC)
               END IF
            END IF
         ELSE
C
C           The proposed working set appears to be linearly dependent.
C
            INFORM = 1
         END IF
      END IF
C
      RETURN
C
C     End of  E04NFR.  (RZADD)
C
99999 FORMAT (/' XXX  Serious ill-conditioning in the working set afte',
     *       'r adding constraint ',I5,/' XXX  Overflow may occur in s',
     *       'ubsequent iterations.',//)
      END
