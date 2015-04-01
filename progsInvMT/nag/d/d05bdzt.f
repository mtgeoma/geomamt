      SUBROUTINE D05BDZ(KC,FC,GC,YS,FORS,NMESH,IORDER,IS,H,TOLNL,WORK,
     *                  LWK,NCT,IENTFN,IFLAG)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     <<<<<<<<<<<<<<<<<<<<<<<<<<<     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C          A code for the solution of nonlinear convolution
C             weakly singular Abel-Volterra equation
C                of the first and the second kind
C
C     --------------------------------------------------------------
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This routine computes the numerical solution of a weakly
C     singular Volterra-Abel integral equation of
C     the form :
C                                  t     -1/2
C     A*y(t) = f(t) + (1/sqrt(pi)) I (t-s)    k(t-s) g(s, y(s)) ds,
C                                  0
C                                                      (t >= 0),
C     where A=0 or A = 1.
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     --------------------------------------------------------------
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H, TOLNL
      INTEGER           IENTFN, IFLAG, IORDER, IS, LWK, NMESH
      CHARACTER         FORS
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(LWK), YS(NMESH+IS)
      INTEGER           NCT(*)
C     .. Function Arguments ..
      DOUBLE PRECISION  FC, GC, KC
      EXTERNAL          FC, GC, KC
C     .. Local Scalars ..
      INTEGER           I1, I10, I2, I3, I4, I5, I6, I7, I8, I9, IFAIL1,
     *                  IP, ISM1, IW, LENP, LFW, LSW, LWKR, N, N1, N12
      CHARACTER         FFT
C     .. External Subroutines ..
      EXTERNAL          D05BDY, D05BYF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, LOG, NINT
C     .. Executable Statements ..
C
      IF (NMESH.LE.32) THEN
         FFT = 'N'
         N = NMESH
      ELSE
         FFT = 'Y'
         N = 32
      END IF
C
      IP = NINT(LOG(DBLE(NMESH))/LOG(DBLE(2)))
      LENP = IP - 1
      N1 = 2**IP
      N12 = 2*N1
      ISM1 = IS - 1
      IW = N12 + 1
      I1 = 1
      I2 = I1 + N - 1
      I3 = I2 + N1 + IS
      I4 = I3 + N12
      I5 = I4 + N1 + 2*ISM1
      I6 = I5 + ISM1*ISM1
      I7 = I6 + ISM1*ISM1
      I8 = I7 + 4*ISM1
      I9 = I8 + NMESH
      I10 = I9 + NMESH
C
      IFAIL1 = 1
      LSW = NMESH + IS
      LFW = 2*NMESH
      LWKR = 2**(LENP+3)
C
      IF (IENTFN.EQ.0) THEN
         CALL D05BYF(IORDER,LENP,LFW,WORK(I8),WORK(I10),LSW,WORK,LWKR,
     *               IFAIL1)
C
         CALL D05BDY(KC,FC,GC,YS(1),N,FFT,LENP,N1,IS,ISM1,H,TOLNL,
     *               WORK(I8),WORK(I10),WORK(I1),WORK(I2),WORK(I3),
     *               WORK(I4),WORK(I9),WORK(I5),WORK(I6),WORK(I7),NCT,
     *               FORS,IFLAG)
C
         IF (IFLAG.NE.0) GO TO 20
C
      ELSE IF (IENTFN.EQ.1) THEN
C
         CALL D05BDY(KC,FC,GC,YS(1),N,FFT,LENP,N1,IS,ISM1,H,TOLNL,
     *               WORK(I8),WORK(I10),WORK(I1),WORK(I2),WORK(I3),
     *               WORK(I4),WORK(I9),WORK(I5),WORK(I6),WORK(I7),NCT,
     *               FORS,IFLAG)
C
         IF (IFLAG.NE.0) GO TO 20
C
      END IF
C
   20 CONTINUE
      RETURN
      END
