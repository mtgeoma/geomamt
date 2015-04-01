      SUBROUTINE D03PRR(XP,UP,IPTS,X,M,U,NPTS,NPDE,IFAIL,C,IXFIX,NFIX)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C***********************************************************************
C  MODIFIED D03PRP FOR UPWIND SCHEME ROUTINES -- CALLS NAG MONOTONICITY
C  PRESERVING INTERPOLATION ROUTINES E01BEF, E01BFF.
C***********************************************************************
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IPTS, M, NFIX, NPDE, NPTS
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NPTS,3), U(NPDE,NPTS), UP(NPDE,IPTS), X(NPTS),
     *                  XP(IPTS)
      INTEGER           IXFIX(*)
C     .. Scalars in Common ..
      INTEGER           IDEV, ITRACE
C     .. Local Scalars ..
      INTEGER           I, IN, IT, J, JJ, JK, JN, K, KK, NKPTS, NXFXP1
C     .. External Subroutines ..
      EXTERNAL          D03PZW, E01BEF, E01BFF
C     .. Common blocks ..
      COMMON            /AD02NM/ITRACE, IDEV
C     .. Save statement ..
      SAVE              /AD02NM/
C     .. Executable Statements ..
      IF (NFIX.EQ.0) THEN
C         No material interfaces present.
         DO 60 J = 1, NPDE
            DO 20 I = 1, NPTS
               C(I,1) = U(J,I)
   20       CONTINUE
            IFAIL = 0
            CALL E01BEF(NPTS,X,C,C(1,2),IFAIL)
            CALL E01BFF(NPTS,X,C,C(1,2),NPTS,XP,C(1,3),IFAIL)
            DO 40 I = 1, NPTS
               UP(J,I) = C(I,3)
   40       CONTINUE
   60    CONTINUE
C
      ELSE
C         Cubic spline routine must be called in piecewise fashion .
         NXFXP1 = NFIX + 1
         JJ = 1
         DO 180 I = 1, NXFXP1
            JK = 0
            IF (I.EQ.1) THEN
               K = 1
            ELSE
               K = IXFIX(I-1) + 0.1D0
            END IF
C            X(K) is the starting point for the next element
            IF (I.EQ.NXFXP1) THEN
               NKPTS = NPTS - K + 1
               KK = NPTS
            ELSE
               NKPTS = IXFIX(I) - K + 1.1D0
               KK = IXFIX(I) + 0.1D0
            END IF
   80       J = JJ + JK
            IF (J.GT.IPTS) THEN
               GO TO 100
            ELSE
               IF (XP(J).LE.X(KK)) THEN
                  JK = JK + 1
C               Counts the no. of new points in this element
                  GO TO 80
               END IF
            END IF
  100       J = J - 1
C           Can interpolate to find new solution values from
C           XP(JJ)th mesh point to XP(J)th mesh point.
C       **************  IF (ITRACE.GE.1)WRITE(IDEV,90)JJ,J,JK,K,NKPTS
C    90          FORMAT(' XP(',I3,') TO XP(',I3,' USING',I4,'POINTS'/
C       ****1********   ' INTERP FROM X(',I3,') USING ',I4,' POINTS')
            IF (JK.GT.0 .AND. NKPTS.GE.4) THEN
C
               DO 160 JN = 1, NPDE
                  DO 120 IN = 1, NPTS
                     C(IN,1) = U(JN,IN)
  120             CONTINUE
                  IFAIL = 0
                  CALL E01BEF(NKPTS,X(K),C(K,1),C(1,2),IFAIL)
                  CALL E01BFF(NKPTS,X(K),C(K,1),C(1,2),NKPTS,XP(JJ),
     *                        C(K,3),IFAIL)
                  DO 140 IN = 1, NKPTS
                     UP(JN,IN+JJ-1) = C(IN+JJ-1,3)
  140             CONTINUE
  160          CONTINUE
C
            ELSE IF (JK.GT.0) THEN
               IT = 1
               CALL D03PZW(XP(JJ),UP(1,JJ),JK,X(K),M,U(1,K),NKPTS,NPDE,
     *                     IT,IFAIL)
            END IF
            JJ = J
  180    CONTINUE
      END IF
      RETURN
      END
