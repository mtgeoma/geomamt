      SUBROUTINE D02PXY(F,NEQ,NWANT,Y,YP,YOLD,YPOLD,STAGES,CALSTG,
     *                  XSTAGE,YTEMP,P)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C$$$$ SUBROUTINE FORMI $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C************************************************
C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
C************************************************
C
C  Purpose:    Forms an interpolating polynomial for use with
C              METHDs 1 or 2.
C
C  Input:      NEQ, NWANT, T, Y(*), YP(*), HOLD, YOLD(*), YPOLD(*),
C              STAGES(NEQ,*), CALSTG
C  Output:     P(*), YTEMP(*), XSTAGE(NEQ)
C  External:   F
C
C  Common:     Initializes:    none
C              Reads:          /DD02PD/ A(*,*), C(*), R(*), METHD, MINTP
C                              /BD02PD/ T, TOLD, HOLD
C              Alters:         /BD02PD/ NFCN
C
C  Comments:
C  =========
C  The integration has reached T with a step HOLD from TOLD = T-HOLD.
C  Y(*),YP(*) and YOLD(*),YPOLD(*) approximate the solution and its
C  derivative at T and TOLD respectively.  STAGES(NEQ,*) holds the
C  stages computed in taking this step. In the case of METHD = 2 it is
C  necessary to compute some more stages in this subroutine. CALSTG
C  indicates whether or not the extra stages need to be computed. A(*,*)
C  and C(*) are used in computing these stages. The extra stages are
C  stored in STAGES(NEQ,*) and XSTAGE(*).  The coefficients of the
C  interpolating polynomials for the first NWANT components of the
C  solution are returned in the array P(*). The polynomial is of degree
C  MINTP = 3 for METHD = 1 and of degree MINTP = 6 for METHD = 2. The
C  vector R(*) is used for workspace when METHD = 2.
C
C     .. Scalar Arguments ..
      INTEGER           NEQ, NWANT
      LOGICAL           CALSTG
C     .. Array Arguments ..
      DOUBLE PRECISION  P(NWANT,*), STAGES(NEQ,*), XSTAGE(*), Y(*),
     *                  YOLD(*), YP(*), YPOLD(*), YTEMP(*)
C     .. Subroutine Arguments ..
      EXTERNAL          F
C     .. Scalars in Common ..
      DOUBLE PRECISION  H, HOLD, T, TOLD
      INTEGER           FLSTP, METHD, MINTP, NFCN, NSTAGE, OKSTP, SVNFCN
      LOGICAL           FIRST, INTP, LAST
C     .. Arrays in Common ..
      DOUBLE PRECISION  A(13,13), B(13), BHAT(13), C(13), E(7), R(11,6)
      INTEGER           PTR(13)
C     .. Local Scalars ..
      DOUBLE PRECISION  D1, D2, D3, D4, HYP, HYPOLD
      INTEGER           I, J, K, L
C     .. Common blocks ..
      COMMON            /BD02PD/T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP,
     *                  FLSTP, FIRST, LAST
      COMMON            /DD02PD/A, B, C, BHAT, R, E, PTR, NSTAGE, METHD,
     *                  MINTP, INTP
C     .. Save statement ..
      SAVE              /BD02PD/, /DD02PD/
C     .. Executable Statements ..
C
      IF (METHD.EQ.1) THEN
C
C  METHD = 1.  Use the cubic Hermite interpolant that is is fully
C  specified by the values and slopes at the two ends of the step.
C
         DO 20 L = 1, NWANT
            D1 = Y(L) - YOLD(L)
            HYP = HOLD*YP(L)
            HYPOLD = HOLD*YPOLD(L)
            D2 = HYP - D1
            D3 = D1 - HYPOLD
            D4 = D2 - D3
            P(L,1) = D2 + D4
            P(L,2) = D4
   20    CONTINUE
C
      ELSE
C
C  METHD = 2.
C
         IF (CALSTG) THEN
C
C  Compute the extra stages needed for interpolation using the facts
C  that
C       1. Stage 1 is YPOLD(*).
C       2. Stage i (i>1) is stored in STAGES(1:NEQ,i).
C       3. This pair is FSAL, i.e. STAGES(1:NEQ,7)=YP(1:NEQ), which
C          frees up STAGES(1:NEQ,7) for use by stage 9.
C       4. XSTAGE(1:NEQ) is used for stage 10.
C       5. The coefficient of stage 2 in the interpolant is always 0, so
C          STAGES(1:NEQ,1) is used for stage 11.
C  The vector P(1:NEQ) is used as workspace for computing the stages.
C
            DO 180 I = 9, 11
               DO 40 L = 1, NEQ
                  YTEMP(L) = A(I,1)*YPOLD(L)
   40          CONTINUE
               DO 140 J = 2, I - 1
                  IF (J.LE.7) THEN
                     DO 60 L = 1, NEQ
                        YTEMP(L) = YTEMP(L) + A(I,J)*STAGES(L,J-1)
   60                CONTINUE
                  ELSE IF (J.EQ.8) THEN
                     DO 80 L = 1, NEQ
                        YTEMP(L) = YTEMP(L) + A(I,J)*YP(L)
   80                CONTINUE
                  ELSE IF (J.EQ.9) THEN
                     DO 100 L = 1, NEQ
                        YTEMP(L) = YTEMP(L) + A(I,J)*STAGES(L,7)
  100                CONTINUE
                  ELSE IF (J.EQ.10) THEN
                     DO 120 L = 1, NEQ
                        YTEMP(L) = YTEMP(L) + A(I,J)*XSTAGE(L)
  120                CONTINUE
                  END IF
  140          CONTINUE
               DO 160 L = 1, NEQ
                  YTEMP(L) = YOLD(L) + HOLD*YTEMP(L)
  160          CONTINUE
               IF (I.EQ.9) THEN
                  CALL F(TOLD+C(I)*HOLD,YTEMP,STAGES(1,7))
                  NFCN = NFCN + 1
               ELSE IF (I.EQ.10) THEN
                  CALL F(TOLD+C(I)*HOLD,YTEMP,XSTAGE)
                  NFCN = NFCN + 1
               ELSE
                  CALL F(TOLD+C(I)*HOLD,YTEMP,STAGES(1,1))
                  NFCN = NFCN + 1
               END IF
  180       CONTINUE
         END IF
C
C  Form the coefficients of the interpolating polynomial in its shifted
C  and scaled form.  The transformation from the form in which the
C  polynomial is derived can be somewhat ill-conditioned.  The terms
C  are grouped so as to minimize the errors of the transformation.
C
C  Coefficient of SIGMA**6
         DO 200 L = 1, NWANT
            P(L,5) = R(5,6)*STAGES(L,4) + ((R(10,6)*XSTAGE(L)+R(8,6)
     *               *YP(L))+(R(7,6)*STAGES(L,6)+R(6,6)*STAGES(L,5))) +
     *               ((R(4,6)*STAGES(L,3)+R(9,6)*STAGES(L,7))+(R(3,6)
     *               *STAGES(L,2)+R(11,6)*STAGES(L,1))+R(1,6)*YPOLD(L))
  200    CONTINUE
C
C  Coefficient of SIGMA**5
         DO 220 L = 1, NWANT
            P(L,4) = (R(10,5)*XSTAGE(L)+R(9,5)*STAGES(L,7)) + ((R(7,5)
     *               *STAGES(L,6)+R(6,5)*STAGES(L,5))+R(5,5)*STAGES(L,4)
     *               ) + ((R(4,5)*STAGES(L,3)+R(8,5)*YP(L))+(R(3,5)
     *               *STAGES(L,2)+R(11,5)*STAGES(L,1))+R(1,5)*YPOLD(L))
  220    CONTINUE
C
C  Coefficient of SIGMA**4
         DO 240 L = 1, NWANT
            P(L,3) = ((R(4,4)*STAGES(L,3)+R(8,4)*YP(L))+(R(7,4)
     *               *STAGES(L,6)+R(6,4)*STAGES(L,5))+R(5,4)*STAGES(L,4)
     *               ) + ((R(10,4)*XSTAGE(L)+R(9,4)*STAGES(L,7))+(R(3,4)
     *               *STAGES(L,2)+R(11,4)*STAGES(L,1))+R(1,4)*YPOLD(L))
  240    CONTINUE
C
C  Coefficient of SIGMA**3
         DO 260 L = 1, NWANT
            P(L,2) = R(5,3)*STAGES(L,4) + R(6,3)*STAGES(L,5) + ((R(3,3)
     *               *STAGES(L,2)+R(9,3)*STAGES(L,7))+(R(10,3)*XSTAGE(L)
     *               +R(8,3)*YP(L))+R(1,3)*YPOLD(L)) + ((R(4,3)
     *               *STAGES(L,3)+R(11,3)*STAGES(L,1))+R(7,3)
     *               *STAGES(L,6))
  260    CONTINUE
C
C  Coefficient of SIGMA**2
C
         DO 280 L = 1, NWANT
            P(L,1) = R(5,2)*STAGES(L,4) + ((R(6,2)*STAGES(L,5)+R(8,2)
     *               *YP(L))+R(1,2)*YPOLD(L)) + ((R(3,2)*STAGES(L,2)
     *               +R(9,2)*STAGES(L,7))+R(10,2)*XSTAGE(L)) + ((R(4,2)
     *               *STAGES(L,3)+R(11,2)*STAGES(L,1))+R(7,2)
     *               *STAGES(L,6))
  280    CONTINUE
C
C  Scale all the coefficients by the step size.
C
         DO 320 K = 1, MINTP - 1
            DO 300 L = 1, NWANT
               P(L,K) = HOLD*P(L,K)
  300       CONTINUE
  320    CONTINUE
C
      END IF
C
      RETURN
      END
