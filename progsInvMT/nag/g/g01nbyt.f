      SUBROUTINE G01NBY(CASE,N,A,LDA,C,LDC,ELA,Q,WK,IRANK,ITEM,IMAX)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C       Checks the existence of the expectation to be calculated
C       using theorems 1-3 of Magnus(1989). Based on routine EXIST
C       by Magnus and Pesaran.
C
C     .. Scalar Arguments ..
      INTEGER           IMAX, IRANK, ITEM, LDA, LDC, N
      CHARACTER         CASE
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,N), C(LDC,*), ELA(*), Q(N,N), WK(N*N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AQMAX, CQMAX, EPS, TEMP
      INTEGER           I, IZERO, J, K
      LOGICAL           NULQAQ, NULQCQ, NULQLA
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, X02AJF, Y90TGF
      EXTERNAL          DDOT, X02AJF, Y90TGF
C     .. External Subroutines ..
      EXTERNAL          DSYMM
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      EPS = X02AJF()
      IF (IRANK.EQ.N) THEN
         ITEM = 1
         IMAX = -1
      ELSE
         IZERO = N - IRANK
C
C        Check to see if A*Q is zero
C
         CALL DSYMM('L','L',N,IZERO,1.0D0,A,LDA,Q,N,0.0D0,WK,N)
         AQMAX = Y90TGF('M','G',N,IZERO,WK,N)
         IF (AQMAX.LT.EPS) THEN
            ITEM = 1
            IMAX = -1
         ELSE
C
C           Check to see if Q'*A*Q is zero
C
            DO 40 I = 1, IZERO
               DO 20 J = 1, I
                  TEMP = DDOT(N,Q(1,I),1,WK((J-1)*N+1),1)
                  IF (TEMP.GT.EPS) THEN
                     NULQAQ = .FALSE.
                     GO TO 60
                  END IF
   20          CONTINUE
   40       CONTINUE
            NULQAQ = .TRUE.
   60       CONTINUE
            IF (CASE.EQ.'R' .OR. CASE.EQ.'r') THEN
               IF (NULQAQ) THEN
                  IMAX = IRANK - 1
                  ITEM = 2
               ELSE
                  IMAX = (IRANK-1)/2
                  ITEM = 3
               END IF
            ELSE IF (CASE.EQ.'L' .OR. CASE.EQ.'l') THEN
C
C              Check to see if Q'*ELA is zero
C
               DO 80 K = 1, IZERO
                  TEMP = DDOT(N,Q(1,K),1,ELA,1)
                  IF (ABS(TEMP).GT.EPS) THEN
                     NULQLA = .FALSE.
                     GO TO 100
                  END IF
   80          CONTINUE
               NULQLA = .TRUE.
  100          CONTINUE
               IF (NULQAQ .AND. NULQLA) THEN
                  IMAX = IRANK
                  ITEM = 2
               ELSE IF (NULQAQ) THEN
                  IMAX = IRANK - 1
                  ITEM = 3
               ELSE IF (NULQLA) THEN
                  IMAX = IRANK/2
                  ITEM = 4
               ELSE
                  IMAX = (IRANK-1)/2
                  ITEM = 5
               END IF
            ELSE IF (CASE.EQ.'Q' .OR. CASE.EQ.'q') THEN
C
C              Check to see if C*Q is zero
C
               CALL DSYMM('L','L',N,IZERO,1.0D0,C,LDC,Q,N,0.0D0,WK,N)
               CQMAX = Y90TGF('M','G',N,IZERO,WK,N)
               IF (CQMAX.LT.EPS) THEN
                  IF (NULQAQ) THEN
                     IMAX = IRANK + 1
                     ITEM = 2
                  ELSE
                     IMAX = (IRANK+1)/2
                     ITEM = 5
                  END IF
               ELSE
C
C                 Check to see if Q'*C*Q is zero
C
                  DO 140 I = 1, IZERO
                     DO 120 J = 1, I
                        TEMP = DDOT(N,Q(1,I),1,WK((J-1)*N+1),1)
                        IF (TEMP.GT.EPS) THEN
                           NULQCQ = .FALSE.
                           GO TO 160
                        END IF
  120                CONTINUE
  140             CONTINUE
                  NULQCQ = .TRUE.
  160             CONTINUE
                  IF (NULQAQ) THEN
                     IF (NULQCQ) THEN
                        IMAX = IRANK
                        ITEM = 3
                     ELSE
                        IMAX = IRANK - 1
                        ITEM = 4
                     END IF
                  ELSE
                     IF (NULQCQ) THEN
                        IMAX = IRANK/2
                        ITEM = 6
                     ELSE
                        IMAX = (IRANK-1)/2
                        ITEM = 7
                     END IF
                  END IF
               END IF
            END IF
         END IF
      END IF
      RETURN
      END
