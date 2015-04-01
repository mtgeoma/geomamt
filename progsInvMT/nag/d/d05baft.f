      SUBROUTINE D05BAF(CK,CG,CF,METHOD,IORDER,ALIM,TLIM,YN,ERREST,NOUT,
     *                  TOL,THRESH,WORK,IWK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 14B REVISED. IER-845 (MAR 1990).
C     MARK 17 REVISED. IER-1634 (JUN 1995).
C     (April 1995). Incorporate D05BWF instead of D05BAW 
C
C     <<<<<<<<<<<<<<<<<<<<<<<<<<<     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C          A code for the solution of convolution Volterra
C                 equation of the second kind
C
C     M.S Derakhshan,   January 1989.
C     ---------------------------------------------------------------
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This subroutine computes the numerical solution of the
C     convolution Volterra integeral equation of the second kind
C
C                           t
C     (1)    y(t) = f(t) +  I k(t - s) g(s,y(s))ds,  (t >= a).
C                           a
C
C     The underlying methods for the numerical solution of (1)
C     are the reducible linear multistep formulae of the Adam's
C     and Backward Differentiation (BD) type.
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     ---------------------------------------------------------------
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D05BAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALIM, THRESH, TLIM, TOL
      INTEGER           IFAIL, IORDER, IWK, NOUT
      CHARACTER*1       METHOD
C     .. Array Arguments ..
      DOUBLE PRECISION  ERREST(NOUT), WORK(IWK), YN(NOUT)
C     .. Function Arguments ..
      DOUBLE PRECISION  CF, CG, CK
      EXTERNAL          CF, CG, CK
C     .. Scalars in Common ..
      DOUBLE PRECISION  GTOL
C     .. Local Scalars ..
      DOUBLE PRECISION  H, HH, HI, UROUND
      INTEGER           I, I1, I2, I3, I4, I5, I6, ILIM1, ILIM2,
     *                  INFO, IOPT, IQ, IRK1, IWT1, L1WT, LSTWT, MINIWK,
     *                  NCOMP, NGRID, NOOFXT, NREC, NSNX, NSTART, NW,
     *                  NW1, NXTPOL, LDSW
C     .. Local Arrays ..
      DOUBLE PRECISION  A(6,0:5), B(0:6), THETA(0:6), V(0:9), WKRK(231),
     *                  WKWGT(717)
      CHARACTER*80      P01REC(5)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          D05BAU, D05BAX, D05BAY, D05BAZ, D05BWF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, INT, LOG, DBLE, SQRT
C     .. Common blocks ..
      COMMON            /AD05BA/GTOL
C     .. Executable Statements ..
C
C     ... Check the inputs. ...
C
      UROUND = SQRT(X02AJF())
      INFO = 1
      IF (METHOD.NE.'A' .AND. METHOD.NE.'a' .AND. METHOD.NE.'B' .AND.
     *    METHOD.NE.'b') THEN
         WRITE (P01REC,FMT=99999) METHOD
         NREC = 1
         GO TO 140
      ELSE IF (IORDER.LT.2 .OR. IORDER.GT.6) THEN
         WRITE (P01REC,FMT=99998) IORDER
         NREC = 2
         GO TO 140
      ELSE IF ((METHOD.EQ.'A' .OR. METHOD.EQ.'a') .AND. IORDER.EQ.2)
     *         THEN
         WRITE (P01REC,FMT=99997)
         NREC = 2
         GO TO 140
      ELSE IF ((METHOD.EQ.'B' .OR. METHOD.EQ.'b') .AND. IORDER.EQ.6)
     *         THEN
         WRITE (P01REC,FMT=99997)
         NREC = 2
         GO TO 140
      ELSE IF (ALIM.LT.0.D0) THEN
         WRITE (P01REC,FMT=99996) ALIM
         NREC = 1
         GO TO 140
      ELSE IF (TLIM.LE.ALIM) THEN
         WRITE (P01REC,FMT=99995) TLIM
         NREC = 2
         GO TO 140
      ELSE IF ((TOL.LT.UROUND) .OR. (TOL.GT.1.D0)) THEN
         WRITE (P01REC,FMT=99994) TOL
         NREC = 2
         GO TO 140
      END IF
C
C     ... Some initialisation. ...
C
      IF (METHOD .EQ. 'A') THEN
         NSTART = IORDER - 2
      ELSE 
         NSTART = IORDER - 1
      END IF
C
      MINIWK = 10*NOUT + 6
      NW = (IWK-6)/5
      H = (TLIM-ALIM)/NOUT
C
C     ... Initializing tolerance for the nonlinear solver ...
C
      GTOL = 1.0D-1*TOL
      IF (TOL.LT.1.0D-6) GTOL = TOL
      IF (NOUT.LE.NSTART) THEN
         INFO = 2
         WRITE (P01REC,FMT=99993) NOUT
         NREC = 3
         GO TO 140
      END IF
C
      IF (IWK.LT.MINIWK) THEN
         INFO = 3
         WRITE (P01REC,FMT=99992) MINIWK, IWK
         NREC = 3
         GO TO 140
      END IF
C
      CALL D05BAY(A,B,THETA,IQ)
      IRK1 = IQ*(NSTART+1) + 1
C
C     ... Setting the parameters of the required method. ...
C
      IF (METHOD .EQ. 'A') THEN
         L1WT = IORDER - 1
         LDSW = L1WT + IORDER - 1
      ELSE
         IF (IORDER.EQ.2) THEN
            L1WT = 45
         ELSE IF (IORDER.EQ.3) THEN
            L1WT = 55
         ELSE IF (IORDER.EQ.4) THEN
            L1WT = 75
         ELSE IF (IORDER.EQ.5) THEN
            L1WT = 115
         END IF
         LDSW = L1WT + IORDER 
      END IF
C
      IWT1 = L1WT + 2
      NW1 = NW + 1
      I1 = NW/2 + 2
      I2 = I1 + NW1
      I3 = I2 + NW1
      I4 = I3 + NW1
      I5 = I4 + NW1
      I6 = I5 + NW/2 + 1
C
      IFW = 0
      CALL D05BWF(METHOD,IORDER,WKWGT(1),L1WT+1,LSTWT,WKWGT(IWT1),LDSW,
     +            NSTART+1,IFW)
C
      DO 20 I = 1, NOUT
         YN(I) = 0.D0
         ERREST(I) = 0.D0
   20 CONTINUE
      DO 40 I = 1, I2 - 1
         WORK(I) = 0.D0
   40 CONTINUE
C
      NXTPOL = 1
      INFO = 0
C
      IF (NSTART.EQ.1) THEN
         CALL D05BAX(CK,CG,CF,ALIM,NXTPOL,TOL,THRESH,H,WKWGT(1),
     *               WKWGT(IWT1),LSTWT,WORK(1),WORK(I1),WORK(I2),
     *               WORK(I3),WORK(I4),WORK(I5),NOUT,NW,INFO)
      ELSE
         CALL D05BAU(CK,CG,CF,ALIM,NSTART,L1WT,LSTWT,NXTPOL,TOL,THRESH,
     *               H,WKWGT(1),WKWGT(IWT1),V,WORK(1),WORK(I1),WORK(I2),
     *               WORK(I3),WORK(I4),WORK(I5),A,B,THETA,IQ,WKRK(1),
     *               WKRK(IRK1),NOUT,NW,INFO)
      END IF
C
      IF (INFO.EQ.4 .OR. INFO.EQ.5) GO TO 60
C
      ILIM1 = 2*NSTART
      ILIM2 = NSTART*NXTPOL
C
      IF (ILIM2.NE.ILIM1) THEN
         CALL D05BAZ(CG,CF,CK,ILIM1,ILIM2,ALIM,NSTART,NXTPOL,L1WT,LSTWT,
     *               TOL,THRESH,H,WKWGT(1),WKWGT(IWT1),WORK(1),WORK(I1),
     *               V,WORK(I2),WORK(I3),WORK(I4),WORK(I5),A,B,THETA,IQ,
     *               WKRK(1),WKRK(IRK1),NOUT,NW,INFO)
C
      END IF
C
   60 CONTINUE
      IF (INFO.EQ.4) THEN
         WRITE (P01REC,FMT=99991)
         NREC = 2
         GO TO 140
      ELSE IF (INFO.EQ.5) THEN
         NCOMP = 5*(NXTPOL*NOUT) + 6
         NXTPOL = NXTPOL/2
         NOOFXT = INT(LOG(DBLE(NXTPOL))/LOG(2.D0)+0.1D0)
         WRITE (P01REC,FMT=99990) NOOFXT, NCOMP
         NREC = 5
         GO TO 140
      END IF
C
C     ... Evaluate values of  the kernel and the forcing term  ...
C     ... for the computation of the remaining solution.       ...
C
      HH = H/NXTPOL
      NSNX = NSTART*NXTPOL
      NGRID = NOUT*NXTPOL
      HI = HH*NSNX
C
      DO 80 I = NSNX + 1, NGRID
         HI = HI + HH
         WORK(I2+I) = CF(ALIM+HI)
         WORK(I3+I) = CK(HI)
   80 CONTINUE
C
      ILIM1 = NSNX
      ILIM2 = NGRID
C
      CALL D05BAZ(CG,CF,CK,ILIM1,ILIM2,ALIM,NSTART,NXTPOL,L1WT,LSTWT,
     *            TOL,THRESH,H,WKWGT(1),WKWGT(IWT1),WORK(1),WORK(I1),V,
     *            WORK(I2),WORK(I3),WORK(I4),WORK(I5),A,B,THETA,IQ,
     *            WKRK(1),WKRK(IRK1),NOUT,NW,INFO)
C
      DO 100 I = 1, NOUT
         YN(I) = WORK(I1+I*NXTPOL)
  100 CONTINUE
C
      DO 120 I = 1, NOUT
         ERREST(I) = ABS(WORK(1+I*NXTPOL/2))
  120 CONTINUE
C
      IF (INFO.EQ.4) THEN
         WRITE (P01REC,FMT=99991)
         NREC = 2
      ELSE IF (INFO.EQ.5) THEN
         INFO = 6
         NCOMP = 5*(NXTPOL*NOUT) + 6
         NXTPOL = NXTPOL/2
         NOOFXT = INT(LOG(DBLE(NXTPOL))/LOG(2.D0)+0.1D0)
         WRITE (P01REC,FMT=99990) NOOFXT, NCOMP
         NREC = 5
      END IF
C
  140 IFAIL = P01ABF(IFAIL,INFO,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, METHOD has not been set to A, a, B or b: ',
     *       'METHOD =',A1)
99998 FORMAT (' ** On entry, IORDER has been set to a value which is',
     *       /' ** either greater than 6 or less than 2:  IORDER = ',
     *       I12)
99997 FORMAT (' ** This method does not exist; Either you have set',/
     *    ' ** METHOD=''A'' and IORDER=2  OR  METHOD=''B'' and IORDER=6'
     *       )
99996 FORMAT (' ** On entry, ALIM has been set to a value less than 0:',
     *       ' ALIM =',1P,D12.5)
99995 FORMAT (' ** On entry, TLIM has been set to a value less than or',
     *       ' equal',/' ** to ALIM: TLIM =',1P,D15.5)
99994 FORMAT (' ** On entry, TOL has been set to a value less than squ',
     *       'are root of',/' ** machine precision or greater then 1.0',
     *       ': TOL =',1P,D12.5)
99993 FORMAT (' ** On entry, NMESH has been set to a value less than',/
     *       ' ** or equal to IORDER-2, when METHOD=''A'', or less than'
     *       ,/' ** or equal to IORDER-1, when METHOD=''B'':  NMESH  = '
     *       ,I12)
99992 FORMAT (' ** LWK, the size of the workspace WORK, is too small to'
     *       ,/' ** initiate the computation; the minimum value is:',
     *       I10,/' ** but LWK has been set to: LWK  =',I12)
99991 FORMAT (' ** The exit is due to an error in the NAG nonlinear so',
     *       'lver.',/' ** The solution is not converging.')
99990 FORMAT (' ** The workspace which has been supplied is too small',
     *       /' ** for  the required accuracy. The number of',/' ** ex',
     *       'trapolations, so far,  is :',I10,/' ** If you require on',
     *       'e more extrapolation extend ',/' ** the size of workspac',
     *       'e to : LWK =',I16)
      END
