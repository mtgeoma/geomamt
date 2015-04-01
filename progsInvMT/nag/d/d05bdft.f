      SUBROUTINE D05BDF(CK,CF,CG,INITWT,IORDER,TLIM,TOLNL,NMESH,YN,WORK,
     *                  LWK,NCT,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     <<<<<<<<<<<<<<<<<<<<<<<<<<<<     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C          A code for the solution of nonlinear convolution
C             weakly singular Abel-Volterra equation
C                      of the second kind
C
C     M.S Derakhshan,
C     Mark 16. Nag Copyright, 1991.
C     --------------------------------------------------------------
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This routine computes the numerical solution of a weakly
C     singular Volterra-Abel integral equation of the 2nd kind of
C     the form :
C                                 t     -1/2
C     y(t) = f(t) + (1/sqrt(pi))  I (t-s)    k(t-s) g(s, y(s)) ds,
C                                 0
C                                                      (t >= 0).
C     The solution  YN(i)  approximates  y(t) at  t = (i-1)*H,
C     i = 1, ..., NMESH. The code is based on
C     the BDF fractional linear multistep method of orders
C     4, 5 and 6 and FFT techniques can be used for the
C     computation of the lag term.
C     For a complete description of the parameters, see the
C     routine specification.
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     --------------------------------------------------------------
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D05BDF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TLIM, TOLNL
      INTEGER           IFAIL, IORDER, LWK, NMESH
      CHARACTER         INITWT
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(LWK), YN(NMESH)
      INTEGER           NCT(NMESH/32+1)
C     .. Function Arguments ..
      DOUBLE PRECISION  CF, CG, CK
      EXTERNAL          CF, CG, CK
C     .. Local Scalars ..
      DOUBLE PRECISION  H, UR
      INTEGER           IENTFN, INFO, IS, LWKCHK, NM, NNMM, NREC
      LOGICAL           ENTRY
      CHARACTER         FORS
C     .. Local Arrays ..
      CHARACTER*80      P01REC(5)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          D05BDZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MOD
C     .. Save statement ..
      SAVE              ENTRY
C     .. Data statements ..
      DATA              ENTRY/.TRUE./
C     .. Executable Statements ..
C
      INFO = 0
      IS = 2*IORDER - 1
      NNMM = NMESH - IS
      H = TLIM/DBLE(NMESH-1)
      UR = 10.D0*X02AJF()
      LWKCHK = (7+IS)*NNMM + (IS-1)*(2*IS+5) + IS*IS
C
C     ... Check input parameters ...
C
      NM = NNMM
   20 CONTINUE
C
      IF (MOD(NM,2).NE.0) THEN
         INFO = 1
         WRITE (P01REC,FMT=99991)
         NREC = 1
         GO TO 40
      ELSE
         IF (NM-2.GT.0) THEN
            NM = NM/2
            GO TO 20
         END IF
      END IF
C
      IF (IORDER.LT.4 .OR. IORDER.GT.6) THEN
         INFO = 1
         WRITE (P01REC,FMT=99999) IORDER
         NREC = 2
         GO TO 40
C
      ELSE IF (TLIM.LE.UR) THEN
         INFO = 1
         WRITE (P01REC,FMT=99997) TLIM
         NREC = 2
         GO TO 40
C
      ELSE IF (NMESH.LT.(2*IORDER+1)) THEN
         INFO = 1
         WRITE (P01REC,FMT=99994) NMESH
         NREC = 2
         GO TO 40
C
      ELSE IF (LWK.LT.LWKCHK) THEN
         INFO = 1
         WRITE (P01REC,FMT=99995) LWK, LWKCHK
         NREC = 3
         GO TO 40
C
      ELSE IF (TOLNL.LE.UR) THEN
         INFO = 1
         WRITE (P01REC,FMT=99990)
         NREC = 1
         GO TO 40
      END IF
C
      IF ((INITWT.EQ.'I') .OR. (INITWT.EQ.'i')) THEN
         IENTFN = 0
      ELSE IF ((INITWT.EQ.'S') .OR. (INITWT.EQ.'s')) THEN
         IENTFN = 1
      ELSE
         INFO = 1
         WRITE (P01REC,FMT=99998) INITWT
         NREC = 2
         GO TO 40
      END IF
C
      IF (ENTRY) THEN
         IF (IENTFN.EQ.1) THEN
            INFO = 1
            WRITE (P01REC,FMT=99996)
            NREC = 2
            GO TO 40
         END IF
         ENTRY = .FALSE.
      END IF
C
      FORS = 'S'
C
      CALL D05BDZ(CK,CF,CG,YN,FORS,NNMM,IORDER,IS,H,TOLNL,WORK,LWK,NCT,
     *            IENTFN,INFO)
C
      IF (INFO.NE.0) THEN
         IF (INFO.EQ.2) THEN
            WRITE (P01REC,FMT=99993)
            NREC = 3
         ELSE IF (INFO.EQ.3) THEN
            WRITE (P01REC,FMT=99992)
            NREC = 3
         END IF
      END IF
C
   40 IFAIL = P01ABF(IFAIL,INFO,SRNAME,NREC,P01REC)
C
      RETURN
C
C
99999 FORMAT (' ** On entry, IORDER is either .gt. 6, or .lt. 4 ',/' *',
     *       '* IORDER = ',I12)
99998 FORMAT (' ** On entry, INITWT is not .eq. to I, i, S or s',/' **',
     *       ' INITWT = ',A1)
99997 FORMAT (' ** On entry, TLIM is .le. 10*eps. ',/' ** TLIM =',1P,
     *       D12.5)
99996 FORMAT (' ** The routine was entered with IENTFN = 1. ')
99995 FORMAT (' ** On entry, LWK  is too small.',/' ** LWK = ',I12,
     *       /' ** Whereas it should be ',I12)
99994 FORMAT (' ** On entry, NMESH is .lt. 2*IORDER+1. ',
     *       /' ** NMESH = ',I12)
99993 FORMAT (' ** An error occured when trying to compute the ',/' **',
     *       ' starting values. For more explanation see  ',/' ** the ',
     *       'Fortran Manual. ')
99992 FORMAT (' ** An error occured when trying to compute the ',/' **',
     *       ' the solution at a specific step. For more',/' ** explan',
     *       'ation see the Fortran Manual. ')
99991 FORMAT (' ** On entry, NMESH does not have the form(2**m+2*IORDE',
     *       'R-1).  ')
99990 FORMAT (' ** On entry, TOLNL is .le. 10*eps. ')
      END
