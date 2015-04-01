      SUBROUTINE D05BWF(METHOD,IORDER,OMEGA,NOMG,LENSW,SW,LDSW,NWT,
     *                  IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     <<<<<<<<<<<<<<<<<<<<<<<<<<<     >>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C          A code for generating the quadrature weights
C          associated with a reducible linear multistep
C                  method for solving Volterra
C                         equations
C
C     Mark 16 Release. NAG Copyright 1991.
C     M.S Derakhshan,  May 1991.
C     ------------------------------------------------------
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C      This subroutine generates the weights associated with
C      the Adams formulae of the orders 3 to 6, and the
C      Backward Differentiation Formulae of the orders
C      2 to 5.
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     ------------------------------------------------------
C
C
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D05BWF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IORDER, LDSW, LENSW, NOMG, NWT
      CHARACTER*1       METHOD
C     .. Array Arguments ..
      DOUBLE PRECISION  OMEGA(0:NOMG-1), SW(LDSW,0:NWT-1)
C     .. Local Scalars ..
      INTEGER           I, IOK, IORDM1, IRLIM1, IRLIM2, J, LIMIT,
     *                  NOMGM1, NREC
C     .. Local Arrays ..
      DOUBLE PRECISION  ADSW(9,0:4), ADWT(0:5), ALFA(0:4)
      CHARACTER*80      P01REC(5)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D05BWM, D05BWN
C     .. Executable Statements ..
C
      NOMGM1 = NOMG - 1
      IOK = 0
      IORDM1 = IORDER - 1
C
      IF (METHOD.EQ.'A') THEN
         LENSW = NOMG + IORDER - 2
      ELSE IF (METHOD.EQ.'B') THEN
         LENSW = NOMG + IORDER - 1
      END IF
C
      IF (METHOD.NE.'A' .AND. METHOD.NE.'a' .AND. METHOD.NE.'B' .AND.
     *    METHOD.NE.'b') THEN
         WRITE (P01REC,FMT=99999) METHOD
         NREC = 1
         IOK = 1
         GO TO 260
C
      ELSE IF (IORDER.LT.2 .OR. IORDER.GT.6) THEN
         WRITE (P01REC,FMT=99998) IORDER
         NREC = 2
         IOK = 2
         GO TO 260
C
      ELSE IF (NOMG.LT.1) THEN
         WRITE (P01REC,FMT=99994) NOMG
         NREC = 1
         IOK = 2
         GO TO 260
C
      ELSE IF ((METHOD.EQ.'A' .OR. METHOD.EQ.'a') .AND. IORDER.EQ.2)
     *         THEN
         WRITE (P01REC,FMT=99997)
         NREC = 2
         IOK = 3
         GO TO 260
C
      ELSE IF ((METHOD.EQ.'B' .OR. METHOD.EQ.'b') .AND. IORDER.EQ.6)
     *         THEN
         WRITE (P01REC,FMT=99997)
         NREC = 2
         IOK = 3
         GO TO 260
C
      ELSE IF ((METHOD.EQ.'A' .OR. METHOD.EQ.'a')
     *         .AND. (NWT.NE.IORDER-1)) THEN
         WRITE (P01REC,FMT=99996)
         NREC = 3
         IOK = 4
         GO TO 260
C
      ELSE IF ((METHOD.EQ.'B' .OR. METHOD.EQ.'b') .AND. (NWT.NE.IORDER))
     *         THEN
         WRITE (P01REC,FMT=99995)
         NREC = 3
         IOK = 4
         GO TO 260
C
      ELSE IF ((METHOD.EQ.'A' .OR. METHOD.EQ.'a')
     *         .AND. (LDSW.LT.(NOMG+IORDER-2))) THEN
         WRITE (P01REC,FMT=99996)
         NREC = 3
         IOK = 5
         GO TO 260
C
      ELSE IF ((METHOD.EQ.'B' .OR. METHOD.EQ.'b')
     *         .AND. (LDSW.LT.(NOMG+IORDER-1))) THEN
         WRITE (P01REC,FMT=99995)
         NREC = 3
         IOK = 5
         GO TO 260
C
      END IF
C
      IF (METHOD.EQ.'A') THEN
C
         CALL D05BWM(IORDER,ADWT,ADSW)
C
         IF (NOMG-1.LT.IORDER-1) THEN
            DO 20 I = 0, NOMG - 1
               OMEGA(I) = ADWT(I)
   20       CONTINUE
            DO 60 J = 0, IORDER - 2
               DO 40 I = 1, LENSW
                  SW(I,J) = ADSW(I,J)
   40          CONTINUE
   60       CONTINUE
         ELSE
C
            DO 80 I = 0, IORDER - 1
               OMEGA(I) = ADWT(I)
   80       CONTINUE
C
            DO 100 I = IORDER, NOMGM1
               OMEGA(I) = 1.D0
  100       CONTINUE
C
            DO 140 J = 0, IORDER - 2
               DO 120 I = 1, 2*IORDER - 3
                  SW(I,J) = ADSW(I,J)
  120          CONTINUE
  140       CONTINUE
            DO 180 J = 0, IORDER - 2
               DO 160 I = (2*IORDER-2), LENSW
                  SW(I,J) = OMEGA(J)
  160          CONTINUE
  180       CONTINUE
         END IF
C
      ELSE IF (METHOD.EQ.'B') THEN
C
         LIMIT = 47 + (IORDER-2)*10 + (IORDER/4)*10 + (IORDER/5)*30
C
         IF (NOMGM1.GT.LIMIT) THEN
C
            CALL D05BWN(IORDER,LIMIT,ALFA,OMEGA,SW,LDSW)
C
            DO 200 I = LIMIT + 1, NOMGM1
               OMEGA(I) = 1.D0
  200       CONTINUE
C
            IRLIM1 = LIMIT + IORDER
            IRLIM2 = LIMIT + IORDER + 1
            DO 240 J = 0, IORDM1
               DO 220 I = IRLIM2, LENSW
                  SW(I,J) = SW(IRLIM1,J)
  220          CONTINUE
  240       CONTINUE
C
         ELSE
C
            CALL D05BWN(IORDER,NOMGM1,ALFA,OMEGA,SW,LDSW)
C
         END IF
C
      END IF
C
  260 IFAIL = P01ABF(IFAIL,IOK,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, METHOD .neq. A, a, B or b: METHOD =',A1)
99998 FORMAT (' ** On entry, IORDER has been set to a value which is',
     *       /' ** either .gt. 6 or .lt. 2:  IORDER = ',I12)
99997 FORMAT (' ** This method does not exist; Either you have set',/
     *    ' ** METHOD=''A'' and IORDER=2  OR  METHOD=''B'' and IORDER=6'
     *       )
99996 FORMAT (' ** On entry, METHOD has been set to A or a,  ',/' ** b',
     *       'ut either  LDSW or NWT have  been ',
     *       /' ** set incorrectly')
99995 FORMAT (' ** On entry, METHOD has been set to B or b, ',/' ** bu',
     *       't either LDSW or NWT been ',/' ** set incorrectly')
99994 FORMAT (' ** On entry, NOMG is .lt. 1, NOMG =',I12)
      END
