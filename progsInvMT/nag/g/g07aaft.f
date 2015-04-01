      SUBROUTINE G07AAF(N,K,CLEVEL,PL,PU,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 15B REVISED. IER-956 (NOV 1991).
C     MARK 17 REVISED. IER-1664 (JUN 1995).
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G07AAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CLEVEL, PL, PU
      INTEGER           IFAIL, K, N
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, C, HDF, P, PHAT, PP, TOL, XK, XN, XN2, Z
      INTEGER           H, IERROR, IFAIL2, IFAULT, L, NREC
C     .. Local Arrays ..
      DOUBLE PRECISION  ZLG(2), ZSM(2)
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  G01CEF, G01FEF, G01FFF, X02AJF
      INTEGER           P01ABF
      EXTERNAL          G01CEF, G01FEF, G01FFF, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          C02AJF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN, DBLE
C     .. Executable Statements ..
C
      IERROR = 0
      NREC = 1
      IF (N.LT.1) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) N
      ELSE IF (K.LT.0) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99998) K
      ELSE IF (N.LT.K) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99997) N, K
      ELSE IF (CLEVEL.LE.0.0D0 .OR. CLEVEL.GE.1.0D0) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99996) CLEVEL
      ELSE
         NREC = 0
         XN = DBLE(N)
         XK = DBLE(K)
         P = 0.5D0 - CLEVEL/2.0D0
         PP = CLEVEL + P
         H = MIN(K,N-K)
         L = MAX(K,N-K)
C
         IF (L.LT.1000000) THEN
            TOL = MAX(50.0D0*X02AJF(),0.5D-12)
            IF (K.EQ.0) THEN
               PL = 0.0D0
            ELSE
               A = XK
               B = XN - XK + 1.0D0
               IFAIL2 = 1
               PL = G01FEF(P,A,B,TOL,IFAIL2)
            END IF
            IF (K.EQ.N) THEN
               PU = 1.0D0
            ELSE
               A = XK + 1.0D0
               B = XN - XK
               IFAIL2 = 1
               PU = G01FEF(PP,A,B,TOL,IFAIL2)
            END IF
         ELSE IF (L.GE.1000000 .AND. H.LE.1000) THEN
            XN2 = 2.0D0*XN
            TOL = MAX(50.0D0*X02AJF(),0.5D-12)
            IF (H.EQ.K) THEN
C
C              Use inverse gamma function
C
               HDF = XK
               IF (K.EQ.0) THEN
                  PL = 0.0D0
               ELSE
                  IFAULT = 1
                  PL = G01FFF(P,HDF,2.0D0,TOL,IFAULT)/XN2
                  IF (IFAULT.EQ.5) THEN
                     IERROR = 2
                     NREC = 2
                     WRITE (P01REC,FMT=99995)
                     PU = 0.0D0
                     GO TO 20
                  END IF
               END IF
C
               HDF = HDF + 1.0D0
               IFAULT = 1
               PU = G01FFF(PP,HDF,2.0D0,TOL,IFAULT)/XN2
               IF (IFAULT.EQ.5) THEN
                  IERROR = 2
                  NREC = 2
                  WRITE (P01REC,FMT=99995)
                  PL = 0.0D0
               END IF
C
            ELSE
               HDF = XN - XK
               IF (K.EQ.N) THEN
                  PU = 1.0D0
               ELSE
                  IFAULT = 1
                  PU = 1.0D0 - G01FFF(P,HDF,2.0D0,TOL,IFAULT)/XN2
                  IF (IFAULT.EQ.5) THEN
                     IERROR = 2
                     NREC = 2
                     WRITE (P01REC,FMT=99995)
                     PU = 0.0D0
                     PL = 0.0D0
                     GO TO 20
                  END IF
               END IF
C
               HDF = HDF + 1.0D0
               IFAULT = 1
               PL = 1.0D0 - G01FFF(PP,HDF,2.0D0,TOL,IFAULT)/XN2
               IF (IFAULT.EQ.5) THEN
                  IERROR = 2
                  NREC = 2
                  WRITE (P01REC,FMT=99995)
                  PL = 0.0D0
                  PU = 0.0D0
               END IF
            END IF
         ELSE
            IFAIL2 = 1
            Z = G01CEF(PP,IFAIL2)
            A = Z*Z/XN + 1.0D0
            PHAT = XK/XN
            B = -Z*Z/XN - 2.0D0*PHAT
            C = PHAT*PHAT
            IFAIL2 = 1
            CALL C02AJF(A,B,C,ZSM,ZLG,IFAIL2)
            PL = MIN(ZSM(1),ZLG(1))
            PU = MAX(ZSM(1),ZLG(1))
         END IF
CRWB - moved following line lower to ensure pl/pu set before accessed
c      END IF
      PL = MAX(0.0D0,PL)
      PU = MIN(1.0D0,PU)
      END IF
C
   20 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, N.lt.1 : N = ',I16)
99998 FORMAT (' ** On entry, K.lt.0 : K = ',I16)
99997 FORMAT (' ** On entry, N.lt.K : N = ',I16,' and K = ',I16)
99996 FORMAT (' ** On entry, CLEVEL.lt.0.0 or CLEVEL.gt.1.0 : CLEVEL = '
     *       ,D13.5)
99995 FORMAT (' ** When using the relationship with the gamma distribu',
     *       'tion the series to',/'    calculate the gamma probabilit',
     *       'ies has failed to converge')
      END
