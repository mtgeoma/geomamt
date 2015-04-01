      SUBROUTINE C06LBF(F,SIGMA0,SIGMA,B,EPSTOL,MMAX,M,ACOEF,ERRVEC,
     *                  IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15B REVISED. IER-946 (NOV 1991).
C
C     C06LBF computes a set of Laguerre expansion coefficients for the
C     inverse Laplace transform of a user-specified analytic function.
C     The value of the inverse Laplace transform at a specified argument
C     can be obtained by a subsequent call to C06LCF.
C
C     C06LBF is derived from the subroutine MODUL1 in the package WEEKS
C     by B.S. Garbow, G. Giunta, J.N. Lyness and A. Murli, Algorithm
C     662: A Fortran software package for the numerical inversion of the
C     Laplace Transform based on Weeks' method, ACM Trans. Math.
C     Software, 14, pp 171-176 (1988).
C
C     INPUT arguments
C
C     F      - procedure - name of function subprogram for the complex
C                     valued Laplace transform to be inverted.  F must
C                     be declared EXTERNAL in the calling program.
C     sigma0 - real - the abscissa of convergence of the Laplace
C                     transform.
C     sigma  - real - the first parameter of the Laguerre expansion.
C                     If sigma is not greater than sigma0, it defaults
C                     and is reset to (sigma0 + 0.7).
C     b      - real - the second parameter of the Laguerre expansion.
C                     If b is less than 2.0*(sigma - sigma0), it
C                     defaults and is reset to 2.5*(sigma - sigma0).
C     epstol - real - the required absolute uniform pseudo accuracy for
C                     the coefficients and inverse Laplace transform
C                     values.
C     mmax   - integer -  an upper limit on the number of coefficients
C                     to be computed.  Note that the maximum number of
C                     Laplace transform evaluations is (mmax/2 + 2).
C
C     OUTPUT arguments
C
C     m      - integer - the number of coefficients actually computed.
C     acoef  - real(mmax) - the array of Laguerre coefficients.
C     errvec - real(8) - an 8-component vector of diagnostic
C                     information.
C          All components are functions of Laguerre coefficients acoef.
C          (1) = Overall estimate of the pseudo-error = (2) + (3) + (4).
C                  Pseudo-error = absolute error / exp(sigma*t) .
C          (2) = Estimate of the discretisation pseudo-error.
C          (3) = Estimate of the truncation pseudo-error.
C          (4) = Estimate of the conditioning pseudo-error, on the basis
C                  of minimal noise levels in function values.
C          (5) = K - Coefficient of the decay function for acoef.
C          (6) = R - Base of the decay function for acoef.
C                  abs(acoef(j+1)) .le. K/R**j for j .ge. m/2 .
C          (7) = ALPHA - Logarithm of the largest acoef.
C          (8) = BETA - Logarithm of the smallest nonzero acoef.
C     ifail  - integer - the output state parameter, takes one of 6
C                        values:
C           0 => Normal termination, estimated error less than epstol.
C           1 => MMAX < 8.
C           2 => Normal termination, but with estimated error bounds
C                slightly larger than epstol.  Note, however, that the
C                actual errors on the final results may be smaller than
C                epstol as bounds independent of t are pessimistic.
C           3 => The round-off level makes it impossible to achieve the
C                required accuracy.
C           4 => The decay rate of the coefficients is too small.
C                It may improve results to increase mmax.
C           5 => The decay rate of the coefficients is too small and the
C                truncation error too large because of round-off error.
C           6 => No error bounds are returned as the behavior of the
C                coefficients does not enable reasonable prediction.
C                Check the value of sigma0.
C                In this case, (errvec(j),j=1,5) are each set to -1.0.
C         NOTE - When ifail is 3, 4, 5 or 6, changing b and sigma
C                may help.  If not, the method should be abandoned.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06LBF')
      DOUBLE PRECISION  ONE, ZERO, BDEF, SIGDEF
      PARAMETER         (ONE=1.0D0,ZERO=0.0D0,BDEF=2.5D0,SIGDEF=0.7D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  B, EPSTOL, SIGMA, SIGMA0
      INTEGER           IFAIL, M, MMAX
C     .. Array Arguments ..
      DOUBLE PRECISION  ACOEF(MMAX), ERRVEC(8)
C     .. Function Arguments ..
      COMPLEX*16        F
      EXTERNAL          F
C     .. Local Scalars ..
      DOUBLE PRECISION  EALFA, EBETA, EPEST, EPMACH, EPREQ, EPSNOI,
     *                  FACT, R1, R2, RCIRC, RCIRCM, REST, RESTM, T1,
     *                  T2, UNLOG, XKEST
      INTEGER           J, JB, JT, NB, NCODE, NREC, NSBOT, NSDIFF,
     *                  NSTATE, NSTOP, NT
C     .. Local Arrays ..
      DOUBLE PRECISION  AB(3), AT(3), PARAMS(3)
      CHARACTER*80      P01REC(4)
C     .. External Functions ..
      COMPLEX*16        C06LBY
      DOUBLE PRECISION  X02AJF, X02AMF
      INTEGER           P01ABF
      EXTERNAL          C06LBY, X02AJF, X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          C06LBZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, EXP, LOG, MAX, MIN, SQRT
C     .. Executable Statements ..
C
C     Set machine-dependent constants
C
      EPMACH = X02AJF()
      UNLOG = -LOG(X02AMF())
C
C     Check mmax
C
      IF (MMAX.LT.8) THEN
         WRITE (P01REC,FMT=99998) MMAX
         NREC = 1
         NSTATE = 1
         GO TO 240
      END IF
C
C     Set parameters to default values and initialize function value
C     in params(3).
C
      IF (SIGMA.LE.SIGMA0) SIGMA = SIGMA0 + SIGDEF
      PARAMS(1) = SIGMA
      IF (B.LT.2*(SIGMA-SIGMA0)) B = BDEF*(SIGMA-SIGMA0)
      PARAMS(2) = B
      PARAMS(3) = ZERO
C
C     Specify the circle radius for coefficient evaluation, the required
C     accuracy, and call the complex differentiation routine C06LBZ.
C
      RCIRC = EXP(-ONE/MAX(MMAX,1024))
C
C     0.367879 = 1/e
C
      EPREQ = 0.367879D0*EPSTOL
      NCODE = 0
      CALL C06LBZ(C06LBY,ZERO,RCIRC,EPREQ,EPMACH,MMAX,NCODE,EPEST,M,
     *            ACOEF,PARAMS,F)
C
C     Unnormalize coefficients and compute parameters for use in C06LCF.
C
      FACT = ONE
      EALFA = ONE
      EBETA = ONE
C
      DO 20 J = 1, M
         ACOEF(J) = ACOEF(J)/FACT
         FACT = RCIRC*FACT
         IF (ACOEF(J).NE.ZERO) THEN
            EALFA = MAX(EALFA,ABS(ACOEF(J)))
            EBETA = MIN(EBETA,ABS(ACOEF(J)))
         END IF
   20 CONTINUE
      ERRVEC(7) = LOG(EALFA)
      ERRVEC(8) = LOG(MAX(EBETA,EPMACH))
C
C     The primary purpose of this routine, that of determining m and
C     computing the coefficients acoef, has now been completed.  The
C     rest of the routine is devoted to calculating the error estimate
C     array errvec and the return code ifail.
C
C     Compute Cauchy inequality constants K=xkest and R=rest in
C     two stages: first from fit to coefficients acoef on their full
C     range, then from fit on the second half only of their range.
C     From these two fits estimate R and then use it to compute K.
C
      NSBOT = 0
      NSTOP = M - 1
C
   40 CONTINUE
C
C     nsbot is lowest s-index and nstop is highest s-index for fit.
C     Index s corresponds to acoef(s+1).
C
      DO 60 J = 1, 3
         IF (ACOEF(NSBOT+J).NE.ZERO) THEN
            AB(J) = MAX(LOG(ABS(ACOEF(NSBOT+J))),-UNLOG)
         ELSE
            AB(J) = -UNLOG
         END IF
         IF (ACOEF(NSTOP-2+J).NE.ZERO) THEN
            AT(J) = MAX(LOG(ABS(ACOEF(NSTOP-2+J))),-UNLOG)
         ELSE
            AT(J) = -UNLOG
         END IF
   60 CONTINUE
      NSDIFF = NSTOP - NSBOT
C
C     Pivot selection - consider nine curves in a natural order.
C
      IF (NSTOP-NSBOT.LT.5) THEN
C
C        Default option when there are not six distinct points.
C
         NB = 1
         NT = 3
         GO TO 180
      END IF
C
      DO 160 NB = 1, 3
         DO 140 NT = 1, 3
            NSDIFF = NSTOP - 2 - NSBOT + NT - NB
C
C           The (nb,nt) curve is fixed at nsbot-1+nb and nstop-3+nt.
C           In the following two loops, the ordinates of this curve at
C           each of the other four points (t1) are compared with the
C           corresponding actual ordinates (t2).  If this curve passes
C           below one of the actual ordinates it is abandoned and the
C           next curve is considered.
C
            DO 80 JB = 1, 3
               IF (JB.NE.NB) THEN
                  T1 = (JB-NB)*AT(NT) + (NSDIFF-(JB-NB))*AB(NB)
                  T2 = NSDIFF*AB(JB)
                  IF (T1.LT.T2) GO TO 120
               END IF
   80       CONTINUE
C
            DO 100 JT = 1, 3
               IF (JT.NE.NT) THEN
                  T1 = (NT-JT)*AB(NB) + (NSDIFF-(NT-JT))*AT(NT)
                  T2 = NSDIFF*AT(JT)
                  IF (T1.LT.T2) GO TO 120
               END IF
  100       CONTINUE
C
            GO TO 180
C
  120       CONTINUE
  140    CONTINUE
  160 CONTINUE
C
C     End of pivot selection.
C
  180 CONTINUE
      R2 = EXP((AB(NB)-AT(NT))/NSDIFF)
C
      IF (NSBOT.EQ.0) THEN
         R1 = R2
         NSBOT = M/2 - 1
         GO TO 40
      END IF
C
C     If both estimates of R are smaller than or too close to one,
C     return with no estimates of K or errors computed.
C
      IF (R1.LT.1.0001D0 .AND. R2.LT.1.0001D0) THEN
         NSTATE = 6
         DO 200 J = 1, 5
            ERRVEC(J) = -ONE
  200    CONTINUE
         ERRVEC(6) = R2
         WRITE (P01REC,FMT=99999) SIGMA0
         NREC = 1
         GO TO 240
      END IF
C
      REST = MAX(R2,(R1+1)/2)
      XKEST = ABS(ACOEF(NSBOT+NB))*REST**(NSBOT-1+NB)
C
      FACT = REST**(NSBOT+3)
      DO 220 J = NSBOT + 3, NSTOP - 3
         XKEST = MAX(XKEST,ABS(ACOEF(J+1))*FACT)
         FACT = REST*FACT
  220 CONTINUE
C
C     Constants K and R are now available.  We may now calculate the
C     rest of the errvec array.
C
      RESTM = REST**M
      ERRVEC(3) = (REST/(REST-1))*XKEST/RESTM
      RCIRCM = RCIRC**M
      ERRVEC(2) = ERRVEC(3)*RCIRCM*(RESTM-1)/(RESTM-RCIRCM)
C
C     Determine eps noise using the average of the
C     (m/2 + 2) absolute values of C06LBY.
C
      EPSNOI = 2*EPMACH*PARAMS(3)/(M+4)
C
C     0.577350 = 1/sqrt(3)
C
      ERRVEC(4) = 0.577350D0*EPSNOI*RCIRC/RCIRCM*SQRT((RCIRCM**2-1)
     *            /(RCIRC**2-1))
      ERRVEC(1) = ERRVEC(2) + ERRVEC(3) + ERRVEC(4)
      ERRVEC(5) = XKEST
      ERRVEC(6) = REST
C
C     Finally, store the output state parameter ifail.
C
      IF (ERRVEC(1).GT.EPSTOL) THEN
         NSTATE = NCODE
         IF (NSTATE.LT.0) NSTATE = 2 - NSTATE
         NSTATE = NSTATE + 1
         IF (NSTATE.EQ.2) THEN
            WRITE (P01REC,FMT=99997)
            NREC = 1
         ELSE IF (NSTATE.EQ.3) THEN
            WRITE (P01REC,FMT=99996)
            NREC = 1
         ELSE IF (NSTATE.EQ.4) THEN
            WRITE (P01REC,FMT=99995) MMAX
            NREC = 2
         ELSE IF (NSTATE.EQ.5) THEN
            WRITE (P01REC,FMT=99994) MMAX
            NREC = 3
         END IF
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99993) ERRVEC(1), EPSTOL
      ELSE
         NSTATE = 0
      END IF
C
  240 CONTINUE
      IFAIL = P01ABF(IFAIL,NSTATE,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** Error bounds cannot be predicted. Check SIGMA0. SIG',
     *       'MA0 =',1P,D14.6)
99998 FORMAT (' ** On entry M .lt. 8: M =',I16)
99997 FORMAT (' ** The estimate of the pseudo-error is slightly larger',
     *       ' than EPSTOL.')
99996 FORMAT (' ** The round-off error level is larger than EPSTOL. In',
     *       'creasing EPSTOL may help.')
99995 FORMAT (' ** The decay rate of the coefficients is too small. In',
     *       'creasing MMAX may help.',/'    MMAX =',I8)
99994 FORMAT (' ** The decay rate of the coefficients is too small and',
     *       /'    round-off error is such that the required accuracy ',
     *       'cannot be obtained.',/'    Increasing MMAX or EPSTOL may',
     *       ' help. MMAX =',I8)
99993 FORMAT ('    Pseudo-error estimate ERRVEC(1) =',1P,D14.6,'  EPST',
     *       'OL =',D14.6)
      END
