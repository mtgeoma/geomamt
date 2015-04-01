      SUBROUTINE D02TKZ(KCOL,RHO,COEF,ACOL,B,ASAVE,M,NEQ,MAXORD,WGTERR,
     *                  JTOL,LTOL,ROOT,WGTMSH,MSTAR,NTOL,TOLIN)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C**********************************************************************
C
C   Purpose
C      Assign (once) values to various array constants.
C
C   Arrays assigned during compilation:
C     CNSTS1 - weights for extrapolation error estimate
C     CNSTS2 - weights for mesh selection
C              (the above weights come from the theoretical form for
C              the collocation error -- see [3])
C
C   Arguments:
C      KCOL   - number of collocation points
C      RHO    - the relative position of the (Gaussian) collocation
C               points on (0,1)
C      COEF   - mmmm!
C      ACOL   - rk-basis coefficients at the collocation points
C      B      -
C      ASAVE  - rk-basis coefficients at points sampled for error
C               estimate
C      M      - orders of the ODEs
C      NEQ    - number of ODEs
C      MAXORD - maximal order of ODEs
C      WGTERR - the particular values of CNSTS1 used for current run
C               (depending on KCOL, M), used in error estimate
C      JTOL   - components of differential system to which tolerances
C               refer (viz, if ltol(i) refers to a derivative of u(j),
C               then jtol(i)=j)
C      LTOL   - components of z(u(x)) for error measure and control
C      ROOT   - reciprocals of expected rates of convergence of compo-
C               nents of z(j) for which tolerances are specified
C      WGTMSH - gotten from the values of CNSTS2 which in turn are
C               the constants in the theoretical expression for the
C               errors; the quantities in WGTMSH are 10 * the values
C               in CNSTS2 so that the mesh selection algorithm
C               is aiming for errors 0.1 * as large as the user
C               requested tolerances.
C      MSTAR  - sum(M(i),i=1,NEQ) number of unknowns in z(u(x))
C      NTOL   - number of tolerances
C      TOLIN  - the tolerances
C
C   Author:
C      R.W. Brankin, NAG Ltd., August 1994
C      (modified version of COLNEW routine CONSTS)
C
C**********************************************************************
C
C
C Assign weights for error estimate
C
C     .. Scalar Arguments ..
      INTEGER           KCOL, MAXORD, MSTAR, NEQ, NTOL
C     .. Array Arguments ..
      DOUBLE PRECISION  ACOL(28,7), ASAVE(28,4), B(28), COEF(KCOL,KCOL),
     *                  RHO(KCOL), ROOT(NTOL), TOLIN(NTOL),
     *                  WGTERR(MSTAR), WGTMSH(NTOL)
      INTEGER           JTOL(NTOL), LTOL(NTOL), M(NEQ)
C     .. Local Scalars ..
      INTEGER           I, IZ, J, JCOMP, KOFF, L, LTOLI, MJ, MTOT
C     .. Local Arrays ..
      DOUBLE PRECISION  CNSTS1(28), CNSTS2(28), DUMMY(1)
C     .. External Subroutines ..
      EXTERNAL          D02TKX, D02TKY
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Data statements ..
      DATA              CNSTS1/.25D0, .625D-1, 7.2169D-2, 1.8342D-2,
     *                  1.9065D-2, 5.8190D-2, 5.4658D-3, 5.3370D-3,
     *                  1.8890D-2, 2.7792D-2, 1.6095D-3, 1.4964D-3,
     *                  7.5938D-3, 5.7573D-3, 1.8342D-2, 4.673D-3,
     *                  4.150D-4, 1.919D-3, 1.468D-3, 6.371D-3,
     *                  4.610D-3, 1.342D-4, 1.138D-4, 4.889D-4,
     *                  4.177D-4, 1.374D-3, 1.654D-3, 2.863D-3/
      DATA              CNSTS2/1.25D-1, 2.604D-3, 8.019D-3, 2.170D-5,
     *                  7.453D-5, 5.208D-4, 9.689D-8, 3.689D-7,
     *                  3.100D-6, 2.451D-5, 2.691D-10, 1.120D-9,
     *                  1.076D-8, 9.405D-8, 1.033D-6, 5.097D-13,
     *                  2.290D-12, 2.446D-11, 2.331D-10, 2.936D-9,
     *                  3.593D-8, 7.001D-16, 3.363D-15, 3.921D-14,
     *                  4.028D-13, 5.646D-12, 7.531D-11, 1.129D-9/
C     .. Executable Statements ..
C
      KOFF = KCOL*(KCOL+1)/2
      IZ = 1
      DO 40 J = 1, NEQ
         MJ = M(J)
         DO 20 L = 1, MJ
            WGTERR(IZ) = CNSTS1(KOFF-MJ+L)
            IZ = IZ + 1
   20    CONTINUE
   40 CONTINUE
C
C Assign array values for mesh selection: wgtmsh, jtol, and root
C
      JCOMP = 1
      MTOT = M(1)
      DO 80 I = 1, NTOL
         LTOLI = LTOL(I)
   60    CONTINUE
         IF (LTOLI.GT.MTOT) THEN
            JCOMP = JCOMP + 1
            MTOT = MTOT + M(JCOMP)
            GO TO 60
         END IF
         JTOL(I) = JCOMP
         WGTMSH(I) = 1.D1*CNSTS2(KOFF+LTOLI-MTOT)/TOLIN(I)
         ROOT(I) = 1.D0/DBLE(KCOL+MTOT-LTOLI+1)
   80 CONTINUE
C
C Specify collocation points
C
      GO TO (100,120,140,160,180,200,220) KCOL
  100 CONTINUE
      RHO(1) = 0.D0
      GO TO 240
  120 CONTINUE
      RHO(2) = .57735026918962576451D0
      RHO(1) = -RHO(2)
      GO TO 240
  140 CONTINUE
      RHO(3) = .77459666924148337704D0
      RHO(2) = .0D0
      RHO(1) = -RHO(3)
      GO TO 240
  160 CONTINUE
      RHO(4) = .86113631159405257523D0
      RHO(3) = .33998104358485626480D0
      RHO(2) = -RHO(3)
      RHO(1) = -RHO(4)
      GO TO 240
  180 CONTINUE
      RHO(5) = .90617984593866399280D0
      RHO(4) = .53846931010568309104D0
      RHO(3) = .0D0
      RHO(2) = -RHO(4)
      RHO(1) = -RHO(5)
      GO TO 240
  200 CONTINUE
      RHO(6) = .93246951420315202781D0
      RHO(5) = .66120938646626451366D0
      RHO(4) = .23861918608319690863D0
      RHO(3) = -RHO(4)
      RHO(2) = -RHO(5)
      RHO(1) = -RHO(6)
      GO TO 240
  220 CONTINUE
      RHO(7) = .949107991234275852452D0
      RHO(6) = .74153118559939443986D0
      RHO(5) = .40584515137739716690D0
      RHO(4) = 0.D0
      RHO(3) = -RHO(5)
      RHO(2) = -RHO(6)
      RHO(1) = -RHO(7)
  240 CONTINUE
C
C...  map (-1,1) to (0,1) by  t = .5 * (1. + x)
C
      DO 260 J = 1, KCOL
         RHO(J) = .5D0*(1.D0+RHO(J))
  260 CONTINUE
C
C...  now find runge-kutta coeffitients b, acol and asave
C...  the values of asave are to be used in  newmsh  and errchk .
C
      DO 300 J = 1, KCOL
         DO 280 I = 1, KCOL
            COEF(I,J) = 0.D0
  280    CONTINUE
         COEF(J,J) = 1.D0
         CALL D02TKX(RHO,COEF(1,J),KCOL)
  300 CONTINUE
      CALL D02TKY(1.D0,COEF,KCOL,MAXORD,B,.FALSE.,DUMMY)
      DO 320 I = 1, KCOL
         CALL D02TKY(RHO(I),COEF,KCOL,MAXORD,ACOL(1,I),.FALSE.,DUMMY)
  320 CONTINUE
      CALL D02TKY(1.D0/6.D0,COEF,KCOL,MAXORD,ASAVE(1,1),.FALSE.,DUMMY)
      CALL D02TKY(1.D0/3.D0,COEF,KCOL,MAXORD,ASAVE(1,2),.FALSE.,DUMMY)
      CALL D02TKY(2.D0/3.D0,COEF,KCOL,MAXORD,ASAVE(1,3),.FALSE.,DUMMY)
      CALL D02TKY(5.D0/6.D0,COEF,KCOL,MAXORD,ASAVE(1,4),.FALSE.,DUMMY)
      RETURN
      END
