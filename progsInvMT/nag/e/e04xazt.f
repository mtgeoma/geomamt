      SUBROUTINE E04XAZ(DEBUG,DONE,FIRST,EPSA,EPSR,FX,INFORM,ITER,ITMAX,
     *                  CDEST,FDEST,SDEST,ERRBND,F1,F2,H,HOPT,HPHI)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     MARK 16 REVISED. IER-1108 (JUL 1993).
C
C     ******************************************************************
C     E04XAZ  implements algorithm  FD, the method described in
C     Gill, P.E., Murray, W., Saunders, M.A., and Wright, M. H.,
C     Computing Forward-Difference Intervals for Numerical Optimization,
C     Siam Journal on Scientific and Statistical Computing, vol. 4,
C     pp. 310-321, June 1983.
C
C     The procedure is based on finding an interval (HPHI) that
C     produces an acceptable estimate of the second derivative, and
C     then using that estimate to compute an interval that should
C     produce a reasonable forward-difference approximation.
C
C     One-sided difference estimates are used to ensure feasibility with
C     respect to an upper or lower bound on X. If X is close to an upper
C     bound, the trial intervals will be negative. The final interval is
C     always positive.
C
C     E04XAZ has been designed to use a reverse communication
C     control structure, i.e., all evaluations of the function occur
C     outside this routine. The calling routine repeatedly calls  E04XAZ
C     after computing the indicated function values.
C
C     E04XAZ  is similar to subroutine FDCORE described in Report
C     SOL 83-6, Documentation of FDCORE and FDCALC, by P.E. Gill,
C     W. Murray,  M.A. Saunders, and M.H. Wright, Department of
C     Operations Research,  Stanford University, Stanford, California
C     94305, June 1983.
C
C     Systems Optimization Laboratory, Stanford University.
C     Based on Fortran 66 Version 2.1 of  FDCORE  written June 1983.
C     Fortran 77 Version written 25-May-1985.
C     This version of  E04XAZ  dated  27-Oct-1992.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  BNDLO, BNDUP
      PARAMETER         (BNDLO=1.0D-3,BNDUP=1.0D-1)
      DOUBLE PRECISION  ZERO, SIXTH, FOURTH
      PARAMETER         (ZERO=0.0D+0,SIXTH=1.6D-1,FOURTH=2.5D-1)
      DOUBLE PRECISION  HALF, TWO
      PARAMETER         (HALF=5.0D-1,TWO=2.0D+0)
      DOUBLE PRECISION  THREE, FOUR, TEN
      PARAMETER         (THREE=3.0D+0,FOUR=4.0D+0,TEN=1.0D+1)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CDEST, EPSA, EPSR, ERRBND, F1, F2, FDEST, FX, H,
     *                  HOPT, HPHI, SDEST
      INTEGER           INFORM, ITER, ITMAX
      LOGICAL           DEBUG, DONE, FIRST
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  AFDMIN, CDSAVE, ERR1, ERR2, FDCERR, FDEST2,
     *                  FDSAVE, HSAVE, OLDCD, OLDH, OLDSD, RHO, SDCERR,
     *                  SDSAVE
      LOGICAL           CE1BIG, CE2BIG, OVERFL, TE2BIG
C     .. Local Arrays ..
      CHARACTER*80      REC(6)
C     .. External Functions ..
      DOUBLE PRECISION  F06BLF
      EXTERNAL          F06BLF
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
C     .. Save statement ..
      SAVE              CDSAVE, FDSAVE, HSAVE, OLDH, RHO, SDSAVE,
     *                  CE1BIG, CE2BIG, TE2BIG
C     .. Executable Statements ..
C
C     ------------------------------------------------------------------
C     Explanation of local variables...
C
C     BNDLO, BNDUP, and RHO control the logic of the routine.
C     BNDLO and BNDUP are the lower and upper bounds that define an
C     acceptable value of the bound on the relative condition error in
C     the second derivative estimate.
C
C     The scalar RHO is the factor by which the interval is multiplied
C     or divided, and also the multiple of the well-scaled interval
C     that is used as the initial trial interval.
C
C     All these values are discussed in the documentation.
C     ------------------------------------------------------------------
C
      ITER = ITER + 1
C
C     Compute the forward-,  backward-,  central-  and second-order
C     difference estimates.
C
      FDEST = F06BLF(F1-FX,H,OVERFL)
      FDEST2 = F06BLF(F2-FX,TWO*H,OVERFL)
C
      OLDCD = CDEST
      CDEST = F06BLF(FOUR*F1-THREE*FX-F2,TWO*H,OVERFL)
C
      OLDSD = SDEST
      SDEST = F06BLF(FX-TWO*F1+F2,H*H,OVERFL)
C
C     Compute  FDCERR  and  SDCERR,  bounds on the relative condition
C     errors in the first and second derivative estimates.
C
      AFDMIN = MIN(ABS(FDEST),ABS(FDEST2))
      FDCERR = F06BLF(EPSA,HALF*ABS(H)*AFDMIN,OVERFL)
      SDCERR = F06BLF(EPSA,FOURTH*ABS(SDEST)*H*H,OVERFL)
C
      IF (DEBUG) THEN
         WRITE (REC,FMT=99999) ITER, FX, H, F1, FDEST, F2, FDEST2,
     *     CDEST, SDEST, FDCERR, SDCERR
         CALL X04BAY(IPRINT,6,REC)
      END IF
C
C     ==================================================================
C     Select the correct case.
C     ==================================================================
      IF (FIRST) THEN
C        ---------------------------------------------------------------
C        First time through.
C        Check whether SDCERR lies in the acceptable range.
C        ------------------------------------------------------------
         FIRST = .FALSE.
         DONE = SDCERR .GE. BNDLO .AND. SDCERR .LE. BNDUP
         TE2BIG = SDCERR .LT. BNDLO
         CE2BIG = SDCERR .GT. BNDUP
         CE1BIG = FDCERR .GT. BNDUP
C
         IF ( .NOT. CE1BIG) THEN
            HSAVE = H
            FDSAVE = FDEST
            CDSAVE = CDEST
            SDSAVE = SDEST
         END IF
C
         RHO = EPSR**(-SIXTH)/FOUR
         IF (TE2BIG) THEN
C
C           The truncation error may be too big  (same as saying
C           SDCERR is too small).  Decrease the trial interval.
C
            RHO = TEN*RHO
            OLDH = H
            H = H/RHO
         ELSE IF (CE2BIG) THEN
C
C           SDCERR is too large.  Increase the trial interval.
C
            OLDH = H
            H = H*RHO
         END IF
      ELSE IF (CE2BIG) THEN
C        ---------------------------------------------------------------
C        During the last iteration,  the trial interval was
C        increased in order to decrease SDCERR.
C        ---------------------------------------------------------------
         IF (CE1BIG .AND. FDCERR.LE.BNDUP) THEN
            CE1BIG = .FALSE.
            HSAVE = H
            FDSAVE = FDEST
            CDSAVE = CDEST
            SDSAVE = SDEST
         END IF
C
C        If SDCERR is small enough, accept H.  Otherwise,
C        increase H again.
C
         DONE = SDCERR .LE. BNDUP
         IF ( .NOT. DONE) THEN
            OLDH = H
            H = H*RHO
         END IF
      ELSE IF (TE2BIG) THEN
C        ---------------------------------------------------------------
C        During the last iteration,  the interval was decreased in order
C        to reduce the truncation error.
C        ---------------------------------------------------------------
         DONE = SDCERR .GT. BNDUP
         IF (DONE) THEN
C
C           SDCERR has jumped from being too small to being too
C           large.  Accept the previous value of H.
C
            H = OLDH
            SDEST = OLDSD
            CDEST = OLDCD
         ELSE
C
C           Test whether FDCERR is sufficiently small.
C
            IF (FDCERR.LE.BNDUP) THEN
               CE1BIG = .FALSE.
               HSAVE = H
               FDSAVE = FDEST
               CDSAVE = CDEST
               SDSAVE = SDEST
            END IF
C
C           Check whether SDCERR is in range.
C
            DONE = SDCERR .GE. BNDLO
C
            IF ( .NOT. DONE) THEN
C
C              SDCERR is still too small, decrease H again.
C
               OLDH = H
               H = H/RHO
            END IF
         END IF
C
      END IF
C
C     ==================================================================
C     We have either finished or have a new estimate of H.
C     ==================================================================
      IF (DONE) THEN
C
C        Sufficiently good second-derivative estimate found.
C        Compute the optimal interval.
C
         HPHI = ABS(H)
         HOPT = TWO*SQRT(EPSA)/SQRT(ABS(SDEST))
C
C        ERR1 is the error bound on the forward-difference estimate
C        with the final value of H.  ERR2 is the difference of FDEST
C        and the central-difference estimate with HPHI.
C
         ERR1 = HOPT*ABS(SDEST)
         ERR2 = ABS(FDEST-CDEST)
         ERRBND = MAX(ERR1,ERR2)
C
C        Set INFORM = 4  if the forward- and central-difference
C        estimates are not close.
C
         INFORM = 0
         IF (ERRBND.GT.HALF*ABS(FDEST)) INFORM = 4
      ELSE
C        ---------------------------------------------------------------
C        Check whether the maximum number of iterations has been
C        exceeded.  If not, exit.
C        ---------------------------------------------------------------
         DONE = ITER .GE. ITMAX
         IF (DONE) THEN
            IF (CE1BIG) THEN
C
C              FDCERR was never small.  Probably a constant function.
C
               INFORM = 1
               HPHI = HOPT
               FDEST = ZERO
               CDEST = ZERO
               SDEST = ZERO
               ERRBND = ZERO
            ELSE IF (CE2BIG) THEN
C
C              FDCERR was small,  but SDCERR was never small.
C              Probably a linear or odd function.
C
               INFORM = 2
               HPHI = ABS(HSAVE)
               HOPT = HPHI
               FDEST = FDSAVE
               CDEST = CDSAVE
               SDEST = ZERO
               ERRBND = TWO*EPSA/HOPT
            ELSE
C
C              The only remaining case occurs when the second
C              derivative is changing too rapidly for an adequate
C              interval to be found (SDCERR remained small even
C              though H was decreased ITMAX times).
C
               INFORM = 3
               HPHI = ABS(HSAVE)
               HOPT = HPHI
               FDEST = FDSAVE
               CDEST = CDSAVE
               SDEST = SDSAVE
               ERRBND = HOPT*ABS(SDEST)/TWO + TWO*EPSA/HOPT
            END IF
         END IF
      END IF
C
      IF (DEBUG) THEN
         WRITE (REC,FMT=99998) CE1BIG, CE2BIG, TE2BIG
         CALL X04BAF(IPRINT,REC(1))
         IF (DONE) THEN
            WRITE (REC,FMT=99997) INFORM, HOPT, ERRBND
            CALL X04BAF(IPRINT,REC(1))
         END IF
      END IF
C
      RETURN
C
C
C     End of  E04XAZ. (CHCORE)
C
99999 FORMAT (/' //E04XAZ//  ITN ',I3,' FX     H',11X,1P,2D16.6,/' //E',
     *       '04XAZ//  F1      FDEST',14X,1P,2D16.6,/' //E04XAZ//  F2 ',
     *       '     FDEST2',13X,1P,2D16.6,/' //E04XAZ//  CDEST   SDEST',
     *       14X,1P,2D16.6,/' //E04XAZ//  FDCERR  SDCERR',13X,1P,2D16.6)
99998 FORMAT (' //E04XAZ//  CE1BIG  CE2BIG  TE2BIG',5X,3L2)
99997 FORMAT (' //E04XAZ//  INFORM  HOPT    ERRBND',I5,1P,2D16.6)
      END
