      SUBROUTINE E04UCL(LSUMRY,UNITQ,N,NCNLN,NFREE,NZ,LDJ1,LDJ2,LDQ,LDR,
     *                  KX,ALFA,GL1,GL2,QPCURV,CJAC1,CJAC2,CJDX1,CJDX2,
     *                  VIOL1,VIOL2,GQ1,GQ2,HPQ,RPQ,QPMUL,R,PENU,Q,W,Y)
C     MARK 17 RE-ISSUE.  NAG COPYRIGHT 1995.
C
C     ==================================================================
C     E04UCL  computes the BFGS update for the approximate Hessian of
C     the Lagrangian.  If the approximate curvature of the Lagrangian
C     function is negative,  a nonnegative penalty vector PENU(i) of
C     minimum two norm is computed such that the approximate curvature
C     of the augmented Lagrangian will be positive. If no finite penalty
C     vector exists,  the BFGS update is performed with the approximate
C     curvature modified to be a small positive value.
C
C     On entry,  GQ1 and GQ2 contain the transformed objective gradients
C     at X1 and X2,  HPQ contains  R'R(pq), the transformed Hessian
C     times the transformed search direction.  The vectors GQ1 and HPQ
C     are not saved.  If the regular BFGS quasi-Newton update could not
C     be performed, the first character of LSUMRY is loaded with 'M'.
C
C     Apr-92: Update always done. Negative CURVL set to TINYCL.
C     Jul-94: Update skipping restored. Seems more stable.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original Fortran 66 version written April 1984.
C     Level 2 BLAS added 12-June-1986.
C     Level 2 matrix routines added 22-Apr-94.
C     This version of E04UCL dated  08-Jul-94.
C     ==================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION  TOLG
      PARAMETER         (TOLG=1.0D-1)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, GL1, GL2, QPCURV
      INTEGER           LDJ1, LDJ2, LDQ, LDR, N, NCNLN, NFREE, NZ
      LOGICAL           UNITQ
      CHARACTER*5       LSUMRY
C     .. Array Arguments ..
      DOUBLE PRECISION  CJAC1(LDJ1,*), CJAC2(LDJ2,*), CJDX1(*),
     *                  CJDX2(*), GQ1(N), GQ2(N), HPQ(N), PENU(*),
     *                  Q(LDQ,*), QPMUL(*), R(LDR,*), RPQ(N), VIOL1(*),
     *                  VIOL2(*), W(N+NCNLN), Y(N)
      INTEGER           KX(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  DRMAX, DRMIN, RCNDBD, RFROBN, RHODMP, RHOMAX,
     *                  RHONRM, SCALE
      LOGICAL           INCRUN
C     .. Local Scalars ..
      DOUBLE PRECISION  BETA, CURVL, ETA, QI, QMAX, QNORM, RTGTP, RTYTS,
     *                  TEST, TINYCL, TRACE1, TRACE2
      INTEGER           I, J
      LOGICAL           GOTPEN, OVERFL, SSBFGS
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, F06BLF
      EXTERNAL          DNRM2, F06BLF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DSCAL, E04NBV, E04NBW,
     *                  F06FBF, F06FCF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Common blocks ..
      COMMON            /DE04UC/RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
      COMMON            /EE04NB/RCNDBD, RFROBN, DRMAX, DRMIN
C     .. Executable Statements ..
C
C     ------------------------------------------------------------------
C     Set CURVL = (GL2 - GL1)'DX,  the approximate curvature of the
C     Lagrangian along DX.  At first, the curvature is not scaled
C     by the steplength alfa.
C     ------------------------------------------------------------------
      CURVL = GL2 - GL1
      TINYCL = QPCURV*TOLG
      SSBFGS = CURVL .LE. ALFA*TINYCL
      GOTPEN = .FALSE.
C     ------------------------------------------------------------------
C     Test if CURVL is sufficiently positive.  If there are no nonlinear
C     constraints,  no update can be performed.
C     ------------------------------------------------------------------
      IF (CURVL.LT.TINYCL) THEN
         LSUMRY(1:1) = 'Modified BFGS'
C
         IF (NCNLN.GT.0) THEN
            CALL F06FBF(NCNLN,ZERO,PENU,1)
            QMAX = ZERO
            DO 20 I = 1, NCNLN
               QI = CJDX2(I)*VIOL2(I) - CJDX1(I)*VIOL1(I)
               QMAX = MAX(QMAX,QI)
               IF (QI.LE.ZERO) THEN
                  W(I) = ZERO
               ELSE
                  W(I) = QI
               END IF
   20       CONTINUE
C
            QNORM = DNRM2(NCNLN,W,1)
C
            TEST = MAX(TINYCL-CURVL,ZERO)
            BETA = F06BLF(QMAX*TEST,QNORM*QNORM,OVERFL)
            IF (BETA.LT.RHOMAX .AND. .NOT. OVERFL) THEN
               GOTPEN = .TRUE.
C              -------------------------------------
C              A modified update has been found.
C              Compute the penalty parameters PENU.
C              -------------------------------------
               LSUMRY(1:1) = ' '
               BETA = TEST/(QNORM*QNORM)
               DO 40 I = 1, NCNLN
                  QI = W(I)
                  PENU(I) = BETA*QI
                  CURVL = CURVL + BETA*QI*QI
   40          CONTINUE
               IF (CURVL.LT.TINYCL) CURVL = TINYCL
            END IF
         END IF
      END IF
C
      IF (CURVL.GE.TINYCL) THEN
C
C        ---------------------------------------------------------------
C        Compute the difference in the Lagrangian gradient.
C        ---------------------------------------------------------------
C        Update GQ1 to include the (augmented) Lagrangian terms.
C
         IF (NCNLN.GT.0) THEN
            IF (GOTPEN) THEN
               CALL DCOPY(NCNLN,VIOL1,1,W,1)
               CALL F06FCF(NCNLN,PENU,1,W,1)
               CALL DAXPY(NCNLN,(-ONE),QPMUL,1,W,1)
            ELSE
               CALL DCOPY(NCNLN,QPMUL,1,W,1)
               CALL DSCAL(NCNLN,(-ONE),W,1)
            END IF
            CALL DGEMV('T',NCNLN,N,ONE,CJAC1,LDJ1,W,1,ZERO,Y,1)
C
            IF (GOTPEN) THEN
               CALL DCOPY(NCNLN,VIOL2,1,W,1)
               CALL F06FCF(NCNLN,PENU,1,W,1)
               CALL DAXPY(NCNLN,(-ONE),QPMUL,1,W,1)
               CALL DSCAL(NCNLN,(-ONE),W,1)
            ELSE
               CALL DCOPY(NCNLN,QPMUL,1,W,1)
            END IF
            CALL DGEMV('T',NCNLN,N,ONE,CJAC2,LDJ2,W,1,ONE,Y,1)
C
            CALL E04NBW(6,N,NZ,NFREE,LDQ,UNITQ,KX,Y,Q,W)
            CALL DAXPY(N,ONE,Y,1,GQ1,1)
         END IF
C
         CALL DCOPY(N,GQ2,1,Y,1)
         CALL DAXPY(N,(-ONE),GQ1,1,Y,1)
         RTGTP = SQRT(QPCURV)
         RTYTS = SQRT(ALFA*CURVL)
C
         IF (SSBFGS) THEN
            ETA = RTYTS/(RTGTP*ALFA)
         ELSE
            ETA = ONE
         END IF
C
         TRACE1 = DNRM2(N,HPQ,1)/RTGTP
         TRACE2 = DNRM2(N,Y,1)/(RTYTS*ETA)
         RFROBN = ETA*SQRT(ABS((RFROBN-TRACE1)*(RFROBN+TRACE1)
     *            +TRACE2**2))
C
C        ===============================================================
C        Update the Cholesky factor of  Q'HQ.
C        ===============================================================
C        Normalize   RPQ and HPQ.
C
         CALL DSCAL(N,(ONE/RTGTP),RPQ,1)
         CALL DSCAL(N,(ONE/RTGTP),HPQ,1)
C
C        Do the self-scaled or regular BFGS update.
C        Form the vector W = gamma * (GQ2 - GQ1) - beta * R'R*PQ,
C        where  gamma = 1/SQRT( CURV ) = 1/SQRT( (GQ2 - GQ1)'SQ )
C
         IF (SSBFGS) THEN
            DO 60 J = 1, N
               CALL DSCAL(J,ETA,R(1,J),1)
               W(J) = Y(J)/RTYTS - HPQ(J)*ETA
   60       CONTINUE
         ELSE
            DO 80 J = 1, N
               W(J) = Y(J)/RTYTS - HPQ(J)
   80       CONTINUE
         END IF
C
C        Perform the update to  R = R + RPQ*W'.
C        RPQ is overwritten. Arrays  GQ1  and  HPQ  are used to store
C        the sines and cosines defined by the plane rotations.
C
         CALL E04NBV(N,0,N,LDR,N,N,R,HPQ,RPQ,W,GQ1,HPQ)
      END IF
C
C     END OF E04UCL. (NPUPDT)
C
      RETURN
      END
