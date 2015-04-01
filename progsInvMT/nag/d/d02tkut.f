      SUBROUTINE D02TKU(XI,N,Z,DMZ,VALSTR,IFIN,KCOL,NEQ,M,MAXORD,MSTAR,
     *                  ASAVE,ERR,ERREST,WGTERR,TOLIN,LTOL,NTOL,ERMX,
     *                  IERMX,IJERMX)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C**********************************************************************
C
C   Purpose:
C      Determine the error estimates and test to see if the
C      error tolerances are satisfied.
C
C   Arguments:
C      XI     - current mesh points
C      N      - number of intervals
C      Z      - approximate solution on mesh XI
C      DMZ    - the array of m(j)-th derivatives of the current
C               solution at each collocation point in each interval
C      VALSTR - values of the previous solution which are needed
C               for the extrapolation-like error estimate; also used
C               to store solution values in current mesh for potential
C               extrapolation test later
C      IFIN   - a 0-1 variable; on return it indicates whether
C               the error tolerances were satisfied
C      KCOL   - the number of collocation points
C      NEQ    - the number of ODEs
C      M      - the orders of the ODEs
C      MAXORD - the maximal order of the ODEs
C      MSTAR  - the total number of unkowns (sum(m(i)),i=1,n)
C      ASAVE  - the rk-basis functions for fixed points used in the
C               error test
C      ERR    - temporary storage used for error estimates
C      ERREST - storage for error estimates
C      WGTERR - weights used in the extrapolation-like error
C               estimate; the array values are assigned in
C               subroutine  CONSTS.
C      TOLIN  - tolerances to use for the error test
C      LTOL   - index of solution components for which the error test
C               is to be perfomed
C      NTOL   - number of solution components to be used in error test
C      ERMX   - maximum error observed
C      IERMX  - index of first subinterval where ERMX occurred
C      IJERMX - index of component corresponding to ERMX
C
C   Author:
C      R.W. Brankin, NAG Ltd., August 1994
C      (modified version of COLNEW routine ERRCHK)
C
C**********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ERMX
      INTEGER           IERMX, IFIN, IJERMX, KCOL, MAXORD, MSTAR, N,
     *                  NEQ, NTOL
C     .. Array Arguments ..
      DOUBLE PRECISION  ASAVE(28,4), DMZ(NEQ,KCOL,N), ERR(MSTAR),
     *                  ERREST(MSTAR), TOLIN(NTOL), VALSTR(MSTAR,4*N),
     *                  WGTERR(MSTAR), XI(N+1), Z(MSTAR,N)
      INTEGER           LTOL(NTOL), M(NEQ)
C     .. Local Scalars ..
      DOUBLE PRECISION  ERFRC, HI, X
      INTEGER           I, IBACK, J, KNEW, KSTORE, L, LTOLJ
C     .. Local Arrays ..
      DOUBLE PRECISION  DUMMY(1)
C     .. External Subroutines ..
      EXTERNAL          D02TKT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
C
      ERMX = 0.0D0
      IERMX = 0
      IJERMX = 0
C
      IFIN = 1
      DO 20 J = 1, MSTAR
         ERREST(J) = 0.D0
   20 CONTINUE
C
C The error estimates are obtained by combining values of the numerical
C solutions for two meshes.
C For each value of IBACK we will consider the two approximations at 2
C points in each of the new subintervals.
C We work backwards through the subintervals so that new values can be
C stored in VALSTR in case they prove to be needed later for an error
C estimate.
C The routine NEWMSH filled in the needed values of the old solution
C in VALSTR.
C
      DO 120 IBACK = 1, N
         I = N + 1 - IBACK
         KNEW = 4*(I-1) + 3
         KSTORE = 2*(I-1) + 2
         HI = XI(I+1) - XI(I)
         X = XI(I) + HI*2.D0/3.D0
         CALL D02TKT(X,VALSTR(1,KNEW),ASAVE(1,3),DUMMY,XI(I),XI(I+1),
     *               Z(1,I),DMZ(1,1,I),KCOL,NEQ,MAXORD,M,MSTAR,.FALSE.,
     *               .FALSE.,DUMMY)
         DO 40 L = 1, MSTAR
            ERR(L) = WGTERR(L)*ABS(VALSTR(L,KNEW)-VALSTR(L,KSTORE))
   40    CONTINUE
         KNEW = 4*(I-1) + 2
         KSTORE = 2*(I-1) + 1
         X = XI(I) + HI/3.D0
         CALL D02TKT(X,VALSTR(1,KNEW),ASAVE(1,2),DUMMY,XI(I),XI(I+1),
     *               Z(1,I),DMZ(1,1,I),KCOL,NEQ,MAXORD,M,MSTAR,.FALSE.,
     *               .FALSE.,DUMMY)
         DO 60 L = 1, MSTAR
            ERR(L) = ERR(L) + WGTERR(L)*ABS(VALSTR(L,KNEW)
     *               -VALSTR(L,KSTORE))
   60    CONTINUE
C
C Find component-wise maximum error estimate
C
         DO 80 L = 1, MSTAR
            ERREST(L) = MAX(ERREST(L),ERR(L))
   80    CONTINUE
C
C Test whether the tolerance requirements are satisfied in the i-th
C interval, but only if already satisifed up to then.
C
C        IF (IFIN.NE.0) THEN
         DO 100 J = 1, NTOL
            LTOLJ = LTOL(J)
            ERFRC = ERR(LTOLJ)/(ABS(Z(LTOLJ,I))+1.D0)
            IF (ERFRC.GT.ERMX) THEN
               ERMX = ERFRC
               IERMX = I
               IJERMX = J
            END IF
            IF (ERFRC.GT.TOLIN(J)) IFIN = 0
C               IF (ERR(LTOLJ).GT.TOLIN(J)*(ABS(Z(LTOLJ,
C     *             I))+1.D0)) IFIN = 0
  100    CONTINUE
C        END IF
  120 CONTINUE
      RETURN
      END
