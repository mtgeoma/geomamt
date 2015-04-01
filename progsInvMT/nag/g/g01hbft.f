      DOUBLE PRECISION FUNCTION G01HBF(TAIL,N,A,B,XMU,SIG,LDSIG,TOL,WK,
     *                                 LWK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Computes the multivariate normal probabilities.
C
C     .. Parameters ..
      INTEGER                          NMAX, LIWK
      PARAMETER                        (NMAX=10,LIWK=250)
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01HBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 TOL
      INTEGER                          IFAIL, LDSIG, LWK, N
      CHARACTER                        TAIL
C     .. Array Arguments ..
      DOUBLE PRECISION                 A(N), B(N), SIG(LDSIG,N),
     *                                 WK(LWK), XMU(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION                 CONST, L1, L2, RHO, SD1, SD2, U1,
     *                                 U2
      INTEGER                          IND
C     .. Arrays in Common ..
      DOUBLE PRECISION                 C(NMAX,NMAX)
C     .. Local Scalars ..
      DOUBLE PRECISION                 ACC, ALPHA, BIG, D, EPSABS, PROB,
     *                                 TEMP
      INTEGER                          I, IERROR, IFAULT, INFO, IWKAJ,
     *                                 J, LWKAJ, MAXPTS, N2, NPTS, NREC
C     .. Local Arrays ..
      DOUBLE PRECISION                 L(NMAX), U(NMAX)
      INTEGER                          IWK(LIWK)
      CHARACTER*80                     P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION                 G01CEF, G01EAF, G01HAF, G01HBY,
     *                                 G01HBZ, X01AAF
      INTEGER                          P01ABF
      EXTERNAL                         G01CEF, G01EAF, G01HAF, G01HBY,
     *                                 G01HBZ, X01AAF, P01ABF
C     .. External Subroutines ..
      EXTERNAL                         D01AJF, D01FCF, F07FDZ
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, LOG, MAX, MIN, DBLE, SQRT
C     .. Common blocks ..
      COMMON                           /AG01HB/C, CONST, L1, L2, U1, U2,
     *                                 SD1, SD2, RHO, IND
C     .. Executable Statements ..
C
      NREC = 1
      IERROR = 1
      PROB = 0.0D0
      IF (N.LT.1) THEN
         WRITE (P01REC,FMT=99999) N
      ELSE IF (N.GT.10) THEN
         WRITE (P01REC,FMT=99992) N
      ELSE IF (LDSIG.LT.N) THEN
         WRITE (P01REC,FMT=99998) LDSIG, N
      ELSE IF (N.GE.3 .AND. LWK.LT.4*N) THEN
         WRITE (P01REC,FMT=99991) LWK
      ELSE IF (N.LE.2 .AND. LWK.LT.1) THEN
         WRITE (P01REC,FMT=99989) LWK
      ELSE
         IF (N.GT.1) THEN
            IF (TOL.LE.0.0D0) THEN
               WRITE (P01REC,FMT=99997) TOL
               GO TO 160
            END IF
         END IF
         IERROR = 0
         IF (N.EQ.1) THEN
            IF (SIG(1,1).LE.0.0D0) THEN
               IERROR = 3
               WRITE (P01REC,FMT=99996)
               GO TO 160
            END IF
            TEMP = 1.0D0/SQRT(SIG(1,1))
            IF (TAIL.EQ.'L' .OR. TAIL.EQ.'l') THEN
               IFAULT = 0
               PROB = G01EAF('L',(B(1)-XMU(1))*TEMP,IFAULT)
            ELSE IF (TAIL.EQ.'U' .OR. TAIL.EQ.'u') THEN
               IFAULT = 0
               PROB = G01EAF('U',(A(1)-XMU(1))*TEMP,IFAULT)
            ELSE IF (TAIL.EQ.'C' .OR. TAIL.EQ.'c') THEN
               IFAULT = 0
               PROB = G01EAF('L',(B(1)-XMU(1))*TEMP,IFAULT) -
     *                G01EAF('L',(A(1)-XMU(1))*TEMP,IFAULT)
               PROB = MAX(PROB,0.0D0)
            ELSE
               IERROR = 1
               WRITE (P01REC,FMT=99994) TAIL
            END IF
            GO TO 160
         END IF
         DO 20 I = 1, N
            IF (SIG(I,I).LE.0.0D0) THEN
               IERROR = 3
               WRITE (P01REC,FMT=99996)
               GO TO 160
            END IF
            C(I,I) = 1.0D0/SQRT(SIG(I,I))
   20    CONTINUE
C
C        set limits
C
         IF (TAIL.EQ.'C' .OR. TAIL.EQ.'c') THEN
            IND = 0
            DO 40 I = 1, N
               IF (B(I).LE.A(I)) THEN
                  IERROR = 2
                  WRITE (P01REC,FMT=99993) I
                  NREC = 2
                  GO TO 160
               ELSE
                  L(I) = C(I,I)*(A(I)-XMU(I))
                  U(I) = C(I,I)*(B(I)-XMU(I))
               END IF
   40       CONTINUE
         ELSE IF (TAIL.EQ.'L' .OR. TAIL.EQ.'l') THEN
            IND = 1
            IFAULT = 1
            BIG = G01CEF(TOL/DBLE(10*N),IFAULT)
            IF (IFAULT.NE.0) BIG = -10.0D0
            DO 60 I = 1, N
               L(I) = BIG
               U(I) = C(I,I)*(B(I)-XMU(I))
   60       CONTINUE
         ELSE IF (TAIL.EQ.'U' .OR. TAIL.EQ.'u') THEN
            IND = -1
            IFAULT = 1
            BIG = G01CEF(TOL/DBLE(10*N),IFAULT)
            IF (IFAULT.NE.0) BIG = -10.0D0
            DO 80 I = 1, N
               L(I) = C(I,I)*(A(I)-XMU(I))
               U(I) = -BIG
   80       CONTINUE
         ELSE
            IERROR = 1
            WRITE (P01REC,FMT=99994) TAIL
            GO TO 160
         END IF
         L1 = L(N-1)
         L2 = L(N)
         U1 = U(N-1)
         U2 = U(N)
         IF (N.EQ.2) THEN
            RHO = SIG(2,1)*C(1,1)*C(2,2)
            IF (ABS(RHO).GT.1.0D0) THEN
               IERROR = 3
               WRITE (P01REC,FMT=99996)
               GO TO 160
            END IF
            IF (IND.EQ.0) THEN
               IFAULT = 0
               PROB = G01HAF(U1,U2,RHO,IFAULT) + G01HAF(L1,L2,RHO,
     *                IFAULT) - G01HAF(L1,U2,RHO,IFAULT) - G01HAF(U1,L2,
     *                RHO,IFAULT)
               PROB = MIN(1.0D0,MAX(PROB,0.0D0))
            ELSE IF (IND.EQ.1) THEN
               IFAULT = 0
               PROB = G01HAF(U1,U2,RHO,IFAULT)
            ELSE
               IFAULT = 0
               PROB = G01HAF(-L1,-L2,RHO,IFAULT)
            END IF
         ELSE
C
C           put correlation matrix in C
C
            DO 120 J = 1, N
               TEMP = C(J,J)
               C(J,J) = 1.0D0
               DO 100 I = J + 1, N
                  C(I,J) = SIG(I,J)*TEMP*C(I,I)
  100          CONTINUE
  120       CONTINUE
C
C           find Choleski decomposition
C
            CALL F07FDZ('L',N,C,NMAX,INFO)
            IF (INFO.NE.0) THEN
               IERROR = 3
               WRITE (P01REC,FMT=99996)
               GO TO 160
            END IF
C
C           compute mean and variance of conditional bivariate
C           normal distribution
C
            SD1 = C(N-1,N-1)
            SD2 = SQRT(C(N,N-1)*C(N,N-1)+C(N,N)*C(N,N))
            RHO = C(N,N-1)/SD2
C
C           compute constant
C
            D = 0.0D0
            DO 140 I = 1, N - 2
               D = D + LOG(C(I,I))
  140       CONTINUE
            CONST = 0.5D0*(LOG(X01AAF(TEMP)*2.0D0)*DBLE(N-2)) + D
C
C           call intergration routine
C
            IF (N.EQ.3) THEN
               IF (LWK.LE.4*LIWK) THEN
                  LWKAJ = LWK
                  IWKAJ = LWK/4
               ELSE
                  LWKAJ = 4*LIWK
                  IWKAJ = LIWK
               END IF
               EPSABS = 0.0D0
               IFAULT = 1
               CALL D01AJF(G01HBZ,L(1),U(1),EPSABS,TOL,PROB,ACC,WK,
     *                     LWKAJ,IWK,IWKAJ,IFAULT)
               IF (IFAULT.EQ.2) THEN
                  IERROR = 5
                  WRITE (P01REC,FMT=99990) TOL
               ELSE IF (IFAULT.NE.0) THEN
                  IERROR = 4
                  WRITE (P01REC,FMT=99995) ACC
               END IF
            ELSE
               NPTS = 0
               N2 = N - 2
               ALPHA = 2**N2 + 2*N2*N2 + 2*N2 + 1
               MAXPTS = (DBLE(LWK)/DBLE(N)-1.0D0)*ALPHA
               IFAULT = 1
               CALL D01FCF(N2,L,U,NPTS,MAXPTS,G01HBY,TOL,ACC,LWK,WK,
     *                     PROB,IFAULT)
               IF (IFAULT.NE.0) THEN
                  IERROR = 4
                  WRITE (P01REC,FMT=99995) ACC
               END IF
            END IF
         END IF
      END IF
  160 G01HBF = PROB
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, N.lt.1 : N = ',I16)
99998 FORMAT (' ** On entry, LDSIG.lt.N : LDSIG = ',I16,' N = ',I16)
99997 FORMAT (' ** On entry, TOL.le.0.0 : TOL = ',D13.5)
99996 FORMAT (' ** On entry, SIG is not positive-definite')
99995 FORMAT (' ** Full accuracy not achieved, relative accuracy = ',
     *       D10.2)
99994 FORMAT (' ** On entry, TAIL is not valid : TAIL = ',A1)
99993 FORMAT (' ** On entry, the ',I16,'th value in B is less than or ',
     *       'equal to',/'    the corresponding value in A')
99992 FORMAT (' ** On entry, N.gt.10 : N = ',I16)
99991 FORMAT (' ** On entry, LWK.lt.4*N : LWK = ',I16)
99990 FORMAT (' ** Accuracy requested by TOL is too strict: TOL = ',
     *       D13.5)
99989 FORMAT (' ** On entry, LWK.lt.1 : LWK = ',I16)
      END
