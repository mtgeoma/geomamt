      SUBROUTINE E04UCN(FEASQP,N,NCLIN,NCNLN,OBJALF,GRDALF,QPCURV,
     *                  ISTATE,CJDX,CMUL,CS,DLAM,RHO,VIOLN,WORK1,WORK2)
C     MARK 13 RE-ISSUE.  NAG COPYRIGHT 1988.
C     MARK 16 REVISED. IER-1083 (JUL 1993).
C     MARK 17 REVISED. IER-1603 (JUN 1995).
C
C     ******************************************************************
C     E04UCN   computes the value and directional derivative of the
C     augmented Lagrangian merit function.  The penalty parameters
C     RHO(j) are boosted if the directional derivative of the resulting
C     augmented Lagrangian function is not sufficiently negative.  If
C     RHO needs to be increased,  the perturbation with minimum two-norm
C     is found that gives a directional derivative equal to  - p'Hp.
C
C     Systems Optimization Laboratory, Stanford University, California.
C     Original version written  27-May-1985.
C     This version of  E04UCN  dated 14-November-1985.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE, TWO
      PARAMETER         (ZERO=0.0D+0,HALF=0.5D+0,ONE=1.0D+0,TWO=2.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  GRDALF, OBJALF, QPCURV
      INTEGER           N, NCLIN, NCNLN
      LOGICAL           FEASQP
C     .. Array Arguments ..
      DOUBLE PRECISION  CJDX(*), CMUL(*), CS(*), DLAM(*), RHO(*),
     *                  VIOLN(*), WORK1(*), WORK2(*)
      INTEGER           ISTATE(*)
C     .. Scalars in Common ..
      DOUBLE PRECISION  RHODMP, RHOMAX, RHONRM, SCALE
      LOGICAL           INCRUN
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION  PTERM, PTERM2, QNORM, RHO1, RHOI, RHOMIN,
     *                  RHONEW, RTMIN, TSCL
      INTEGER           I, NPLIN
      LOGICAL           BOOST, OVERFL
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2, F06BLF
      EXTERNAL          DDOT, DNRM2, F06BLF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, F06FCF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Common blocks ..
      COMMON            /AX02ZA/WMACH
      COMMON            /DE04UC/RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
C     .. Save statement ..
      SAVE              /AX02ZA/
C     .. Executable Statements ..
C
      IF (NCNLN.EQ.0) RETURN
C
      RTMIN = WMACH(6)
C
      OBJALF = OBJALF - DDOT(NCNLN,CMUL,1,CS,1)
      GRDALF = GRDALF - DDOT(NCNLN,DLAM,1,CS,1)
C
      CALL DCOPY(NCNLN,CS,1,WORK1,1)
C
      IF ( .NOT. FEASQP) THEN
         NPLIN = N + NCLIN
C
         DO 20 I = 1, NCNLN
            IF (ISTATE(NPLIN+I).LT.0 .OR. VIOLN(I).NE.ZERO) WORK1(I)
     *          = -CJDX(I)
   20    CONTINUE
      END IF
C
      GRDALF = GRDALF + DDOT(NCNLN,WORK1,1,CMUL,1)
C
      IF (FEASQP) THEN
C
C        Find the quantities that define  rhomin, the vector of minimum
C        two-norm such that the directional derivative is one half of
C        approximate curvature   - (dx)'H(dx).
C
         DO 40 I = 1, NCNLN
            IF (ABS(CS(I)).LE.RTMIN) THEN
               WORK2(I) = ZERO
            ELSE
               WORK2(I) = CS(I)**2
            END IF
   40    CONTINUE
C
         QNORM = DNRM2(NCNLN,WORK2,1)
         TSCL = F06BLF(GRDALF+HALF*QPCURV,QNORM,OVERFL)
         IF (ABS(TSCL).LE.RHOMAX .AND. .NOT. OVERFL) THEN
C           ------------------------------------------------------------
C           Bounded  RHOMIN  found.  The final value of  RHO(J)  will
C           never be less than  RHOMIN(j).  If the  QP  was feasible,  a
C           trial value  RHONEW  is computed that is equal to the
C           geometric mean of the previous  RHO  and a damped value of
C           RHOMIN.  The new  RHO  is defined as  RHONEW  if it is less
C           than half the previous  RHO  and greater than  RHOMIN.
C           ------------------------------------------------------------
            SCALE = ONE
            DO 60 I = 1, NCNLN
               RHOMIN = MAX((WORK2(I)/QNORM)*TSCL,ZERO)
               RHOI = RHO(I)
C
               RHONEW = SQRT(RHOI*(RHODMP+RHOMIN))
               IF (RHONEW.LT.HALF*RHOI) RHOI = RHONEW
               IF (RHOI.LT.RHOMIN) RHOI = RHOMIN
               RHO(I) = RHOI
   60       CONTINUE
C
            RHO1 = RHONRM
            RHONRM = DNRM2(NCNLN,RHO,1)
C
C           ------------------------------------------------------------
C           If  INCRUN = .TRUE.,  there has been a run of iterations in
C           which the norm of  RHO  has not decreased.  Conversely,
C           INCRUN = false  implies that there has been a run of
C           iterations in which the norm of RHO has not increased.  If
C           INCRUN changes during this iteration the damping parameter
C           RHODMP is increased by a factor of two.  This ensures that
C           RHO(j) will oscillate only a finite number of times.
C           ------------------------------------------------------------
            BOOST = .FALSE.
            IF (INCRUN .AND. RHONRM.LT.RHO1) BOOST = .TRUE.
            IF ( .NOT. INCRUN .AND. RHONRM.GT.RHO1) BOOST = .TRUE.
            IF (BOOST) THEN
               RHODMP = TWO*RHODMP
               INCRUN = .NOT. INCRUN
            END IF
         END IF
C
      ELSE
C
C        The  QP  was infeasible.  Do not alter the penalty parameters,
C        but compute the scale factor so that the constraint violations
C        are reduced.
C
         CALL F06FCF(NCNLN,RHO,1,WORK1,1)
         PTERM2 = DDOT(NCNLN,WORK1,1,CS,1)
C
         SCALE = RHOMAX
         TSCL = F06BLF(GRDALF,PTERM2,OVERFL)
         IF (TSCL.GT.SCALE .AND. TSCL.LE.RHOMAX/(ONE+RHONRM)
     *       .AND. .NOT. OVERFL) SCALE = TSCL
C
         CALL DCOPY(NCNLN,CS,1,WORK1,1)
      END IF
C
C     ------------------------------------------------------------------
C     Compute the new value and directional derivative of the
C     merit function.
C     ------------------------------------------------------------------
      CALL F06FCF(NCNLN,RHO,1,WORK1,1)
C
      PTERM = DDOT(NCNLN,WORK1,1,CS,1)
      OBJALF = OBJALF + HALF*SCALE*PTERM
C
      IF (FEASQP) PTERM2 = PTERM
C
      GRDALF = GRDALF - SCALE*PTERM2
C
      RETURN
C
C
C     End of  E04UCN. (NPMRT)
C
      END
