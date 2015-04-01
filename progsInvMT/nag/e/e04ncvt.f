      SUBROUTINE E04NCV(UNITQ,INFORM,IFIX,IADD,JADD,NACTIV,NZ,NFREE,
     *                  NRANK,NRES,NGQ,N,LDA,LDZY,LDR,LDT,KX,CONDMX,A,R,
     *                  T,RES,GQ,ZY,W,C,S,MSGLVL)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1072 (JUL 1993).
C     MARK 17 REVISED. IER-1585 (JUN 1995).
C
C     ******************************************************************
C     E04NCV  updates the factorization,  A(free) * (Z Y) = (0 T),  when
C     a constraint is added to the working set.  If  NRANK .gt. 0, the
C     factorization  ( R ) = PCQ  is also updated,  where  C  is the
C                    ( 0 )
C     least squares matrix,  R  is upper-triangular,  and  P  is an
C     orthogonal matrix.  The matrices  C  and  P  are not stored.
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
C     If  NRES .GT. 0,  the row transformations are applied to the rows
C     of the  (N by NRES)  matrix  RES.
C     If  NGQ .GT. 0,  the column transformations are applied to the
C     columns of the  (NGQ by N)  matrix  GQ'.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October--1984.
C     Level-2 matrix routines  added 25-Apr-1988.
C     This version of  E04NCV  dated 28-May-1988.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CONDMX
      INTEGER           IADD, IFIX, INFORM, JADD, LDA, LDR, LDT, LDZY,
     *                  MSGLVL, N, NACTIV, NFREE, NGQ, NRANK, NRES, NZ
      LOGICAL           UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(N), GQ(N,*), R(LDR,*), RES(N,*),
     *                  S(N), T(LDT,*), W(N), ZY(LDZY,*)
      INTEGER           KX(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, DTMAX, DTMIN, EPSPT3, EPSPT5, EPSPT8,
     *                  EPSPT9
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  COND, CONDBD, DTNEW, TDTMAX, TDTMIN
      INTEGER           I, NANEW, NFMIN, NPIV, NT
      LOGICAL           BOUND, OVERFL
C     .. Local Arrays ..
      CHARACTER*80      REC(5)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, F06BLF
      EXTERNAL          DNRM2, F06BLF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DSCAL, E04NBW, F06FLF, F06FQF, F06QHF,
     *                  F06QKF, F06QNZ, F06QRF, F06QVF, F06QXF, F06QZZ,
     *                  X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
C     .. Executable Statements ..
C
C     If the condition estimator of the updated factors is greater than
C     CONDBD,  a warning message is printed.
C
      CONDBD = ONE/EPSPT9
C
      OVERFL = .FALSE.
      BOUND = JADD .LE. N
C
      IF (BOUND) THEN
C        ===============================================================
C        A simple bound has entered the working set.  IADD  is not used.
C        ===============================================================
         NANEW = NACTIV
C
         IF (UNITQ) THEN
C
C           Q  is not stored, but KX defines an ordering of the columns
C           of the identity matrix that implicitly define  Q.
C           Define the sequence of pairwise interchanges P that moves
C           the newly-fixed variable to position NFREE.
C           Reorder KX accordingly.
C
            DO 20 I = 1, NFREE - 1
               IF (I.GE.IFIX) THEN
                  W(I) = I + 1
                  KX(I) = KX(I+1)
               ELSE
                  W(I) = I
               END IF
   20       CONTINUE
C
         ELSE
C           ------------------------------------------------------------
C           Q  is stored explicitly.
C           ------------------------------------------------------------
C           Set  W = the  (IFIX)-th  row of  Q.
C           Move the  (NFREE)-th  row of  Q  to position  IFIX.
C
            CALL DCOPY(NFREE,ZY(IFIX,1),LDZY,W,1)
            IF (IFIX.LT.NFREE) THEN
               CALL DCOPY(NFREE,ZY(NFREE,1),LDZY,ZY(IFIX,1),LDZY)
               KX(IFIX) = KX(NFREE)
            END IF
         END IF
         KX(NFREE) = JADD
      ELSE
C        ===============================================================
C        A general constraint has entered the working set.
C        IFIX  is not used.
C        ===============================================================
         NANEW = NACTIV + 1
C
C        Transform the incoming row of  A  by  Q'.  Use C as workspace.
C
         CALL DCOPY(N,A(IADD,1),LDA,W,1)
         CALL E04NBW(8,N,NZ,NFREE,LDZY,UNITQ,KX,W,ZY,C)
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
C           set. Update the estimate of the condition number.
C
            TDTMAX = MAX(DTNEW,DTMAX)
            TDTMIN = MIN(DTNEW,DTMIN)
            COND = F06BLF(TDTMAX,TDTMIN,OVERFL)
         END IF
C
         IF (COND.GT.CONDMX .OR. OVERFL) GO TO 60
C
         IF (UNITQ) THEN
C
C           First general constraint added.  Set  Q = I.
C
            CALL F06QHF('General',NFREE,NFREE,ZERO,ONE,ZY,LDZY)
            UNITQ = .FALSE.
         END IF
      END IF
C
      IF (BOUND) THEN
         NPIV = NFREE
      ELSE
         NPIV = NZ
      END IF
C
      NT = MIN(NRANK,NPIV)
C
      IF (UNITQ) THEN
C        ---------------------------------------------------------------
C        Q (i.e., ZY) is not stored explicitly.
C        Apply the sequence of pairwise interchanges P that moves the
C        newly-fixed variable to position NFREE.
C        ---------------------------------------------------------------
         IF (NGQ.GT.0) CALL F06QKF('Left','Transpose',NFREE-1,W,NGQ,GQ,
     *                             N)
C
         IF (NRANK.GT.0) THEN
C
C           Apply the pairwise interchanges to the triangular part of R.
C           The subdiagonal elements generated by this process are
C           stored in  s(1), s(2), ..., s(nt-1).
C
            CALL F06QNZ('Right',N,IFIX,NT,S,R,LDR)
C
            IF (NT.LT.NPIV) THEN
C
C              R is upper trapezoidal.  Apply the interchanges in
C              columns  nt  thru  npiv.
C
               DO 40 I = IFIX, NT - 1
                  W(I) = I
   40          CONTINUE
C
               CALL F06QKF('Right','Normal',NFREE-1,W,NT,R,LDR)
            END IF
C
C           Eliminate the subdiagonal elements of R with a left-hand
C           sweep of rotations P2 in planes (1,2), (2,3), ...,(nt-1,nt).
C           Apply P2 to RES.
C
            CALL F06QRF('Left ',N,IFIX,NT,C,S,R,LDR)
            IF (NRES.GT.0) CALL F06QXF('Left','Variable','Forwards',NT,
     *                                 NRES,IFIX,NT,C,S,RES,N)
         END IF
      ELSE
C        ---------------------------------------------------------------
C        Full matrix Q.  Define a sweep of plane rotations P such that
C                           Pw = beta*e(npiv).
C        The rotations are applied in the planes (1,2), (2,3), ...,
C        (npiv-1,npiv).  The rotations must be applied to ZY, R, T
C        and GQ'.
C        ---------------------------------------------------------------
         CALL F06FQF('Varble','Forwrds',NPIV-1,W(NPIV),W,1,C,S)
C
         IF (BOUND .AND. NACTIV.GT.0) THEN
C
            CALL DCOPY(NACTIV,S(NZ),1,W(NZ),1)
C
            S(NZ) = S(NZ)*T(NACTIV,NZ+1)
            T(NACTIV,NZ+1) = C(NZ)*T(NACTIV,NZ+1)
C
            CALL F06QZZ('Create',NACTIV,1,NACTIV,C(NZ+1),S(NZ+1),
     *                  T(1,NZ+1),LDT)
            CALL DCOPY(NACTIV,S(NZ),1,T(NACTIV,NZ),LDT-1)
C
            CALL DCOPY(NACTIV,W(NZ),1,S(NZ),1)
         END IF
C
         IF (NGQ.GT.0) CALL F06QXF('Left ','Variable','Forwards',NPIV,
     *                             NGQ,1,NPIV,C,S,GQ,N)
         CALL F06QXF('Right','Variable','Forwards',NFREE,NFREE,1,NPIV,C,
     *               S,ZY,LDZY)
C
         IF (NRANK.GT.0) THEN
C
C           Apply the rotations to the triangular part of R.
C           The subdiagonal elements generated by this process are
C           stored in  s(1),  s(2), ..., s(nt-1).
C
            NT = MIN(NRANK,NPIV)
            CALL F06QVF('Right',N,1,NT,C,S,R,LDR)
C
            IF (NT.LT.NPIV) THEN
C
C              R is upper trapezoidal.  Pretend R is (nt x n) and
C              apply the rotations in columns  nt  thru  npiv.
C
               CALL F06QXF('Right','Variable','Forwards',NT,N,NT,NPIV,C,
     *                     S,R,LDR)
            END IF
C
C           Eliminate the subdiagonal elements of R with a left-hand
C           sweep of rotations P2 in planes (1,2), (2,3), ...,(nt-1,nt).
C           Apply P2 to RES.
C
            CALL F06QRF('Left ',N,1,NT,C,S,R,LDR)
            IF (NRES.GT.0) CALL F06QXF('Left','Variable','Forwards',NT,
     *                                 NRES,1,NT,C,S,RES,N)
         END IF
C
         IF (BOUND) THEN
C
C           The last row and column of ZY has been transformed to plus
C           or minus the unit vector E(NFREE).  We can reconstitute the
C           columns of GQ and R corresponding to the new fixed variable.
C
            IF (W(NFREE).LT.ZERO) THEN
               NFMIN = MIN(NRANK,NFREE)
               IF (NFMIN.GT.0) CALL DSCAL(NFMIN,-ONE,R(1,NFREE),1)
               IF (NGQ.GT.0) CALL DSCAL(NGQ,-ONE,GQ(NFREE,1),N)
            END IF
C
C           ------------------------------------------------------------
C           The diagonals of T have been altered.  Recompute the
C           largest and smallest values.
C           ------------------------------------------------------------
            IF (NACTIV.GT.0) THEN
               CALL F06FLF(NACTIV,T(NACTIV,NZ),LDT-1,TDTMAX,TDTMIN)
               COND = F06BLF(TDTMAX,TDTMIN,OVERFL)
            END IF
         ELSE
C           ------------------------------------------------------------
C           General constraint.  Install the new row of T.
C           ------------------------------------------------------------
            CALL DCOPY(NANEW,W(NZ),1,T(NANEW,NZ),LDT)
         END IF
      END IF
C
C     ==================================================================
C     Prepare to exit.  Check the magnitude of the condition estimator.
C     ==================================================================
   60 IF (NANEW.GT.0) THEN
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
C
C     End of  E04NCV. (LSADD)
C
99999 FORMAT (/' XXX  Serious ill-conditioning in the working set afte',
     *       'r adding constraint ',I5,/' XXX  Overflow may occur in s',
     *       'ubsequent iterations.',//)
      END
