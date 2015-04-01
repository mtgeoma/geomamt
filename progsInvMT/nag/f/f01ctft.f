      SUBROUTINE F01CTF(TRANSA,TRANSB,M,N,ALPHA,A,LDA,BETA,B,LDB,C,LDC,
     *                  IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 14C REVISED. IER-883 (NOV 1990).
C
C  Purpose
C  =======
C
C  F01CTF performs one of the matrix-matrix operations
C
C     C := alpha*op( A ) + beta*op( B ),
C
C  where  op( X ) is one of
C
C     op( X ) = X   or   op( X ) = X',
C
C  alpha and beta are scalars, and A, B and C are matrices, with op( A )
C  an m by n matrix,  op( B )  an m by n matrix, and C an m by n matrix.
C
C  Parameters
C  ==========
C
C  TRANSA - CHARACTER*1.
C           On entry, TRANSA specifies the form of op( A ) to be used in
C           the matrix addition as follows:
C
C              TRANSA = 'N' or 'n',  op( A ) = A.
C
C              TRANSA = 'T' or 't',  op( A ) = A'.
C
C              TRANSA = 'C' or 'c',  op( A ) = A'.
C
C           Unchanged on exit.
C
C  TRANSB - CHARACTER*1.
C           On entry, TRANSB specifies the form of op( B ) to be used in
C           the matrix addition as follows:
C
C              TRANSB = 'N' or 'n',  op( B ) = B.
C
C              TRANSB = 'T' or 't',  op( B ) = B'.
C
C              TRANSB = 'C' or 'c',  op( B ) = B'.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry,  M  specifies  the number  of rows of the matrices
C           op( A ) and op( B ). M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N  specifies the number of columns of the matrices
C           op( A ) and op( B ). N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - REAL            .
C           On entry, ALPHA specifies the scalar alpha.  When  ALPHA  is
C           supplied as zero then A need not be set on entry.
C           Unchanged on exit.
C
C  A      - REAL             array of DIMENSION ( LDA, ma ), where ma is
C           m when TRANSA = 'N' or 'n', and is n otherwise.
C           If alpha is zero, then A need not be set on entry.
C           If alpha is not zero, then before entry with TRANSA = 'N' or
C           'n', the leading m by n part of the array A must contain the
C           matrix  A, otherwise the leading n by m part of the array  A
C           must contain the matrix A.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
C           LDA must be at least  max( 1, m ), otherwise  LDA must be at
C           least  max( 1, n ).
C           Unchanged on exit.
C
C  BETA   - REAL            .
C           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
C           supplied as zero then B need not be set on entry.
C           Unchanged on exit.
C
C  B      - REAL             array of DIMENSION ( LDB, mb ), where mb is
C           m when TRANSB = 'N' or 'n', and is n otherwise.
C           If beta is zero, then B need not be set on entry.
C           If beta is not zero, then before entry with  TRANSB = 'N' or
C           'n', the leading m by n part of the array B must contain the
C           matrix  B, otherwise the leading n by m part of the array  B
C           must contain the matrix B.
C           Unchanged on exit.
C
C  LDB    - INTEGER.
C           On entry, LDB specifies the first dimension of B as declared
C           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
C           LDB must be at least  max( 1, m ), otherwise  LDB must be at
C           least  max( 1, n ).
C           Unchanged on exit.
C
C  C      - REAL             array of DIMENSION ( LDC, N ).
C           On exit, the array  C  is overwritten by the  m by n  matrix
C           alpha*op( A ) + beta*op( B ).
C
C  LDC    - INTEGER.
C           On entry, LDC specifies the first dimension of C as declared
C           in  the  calling  (sub)  program.   LDC  must  be  at  least
C           max( 1, m ).
C           Unchanged on exit.
C
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01CTF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA, BETA
      INTEGER           IFAIL, LDA, LDB, LDC, M, N
      CHARACTER*1       TRANSA, TRANSB
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,N)
C     .. Local Scalars ..
      INTEGER           I, IERR, J, NREC, NROWA, NROWB
      LOGICAL           NOTA, NOTB
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
C     transposed.
C
      NOTA = TRANSA .EQ. 'N' .OR. TRANSA .EQ. 'n'
      NOTB = TRANSB .EQ. 'N' .OR. TRANSB .EQ. 'n'
      IF (NOTA) THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      IF (NOTB) THEN
         NROWB = M
      ELSE
         NROWB = N
      END IF
C
C     Test the input parameters.
C
      IERR = 0
      IF ( .NOT. NOTA .AND. .NOT. (TRANSA.EQ.'C' .OR. TRANSA.EQ.'c')
     *     .AND. .NOT. (TRANSA.EQ.'T' .OR. TRANSA.EQ.'t')) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99999) 'TRANSA', 'TRANSA', TRANSA
      ELSE IF ( .NOT. NOTB .AND. .NOT. (TRANSB.EQ.'C' .OR. TRANSB.EQ.
     *         'c') .AND. .NOT. (TRANSB.EQ.'T' .OR. TRANSB.EQ.'t')) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99999) 'TRANSB', 'TRANSB', TRANSB
      ELSE IF (M.LT.0 .OR. N.LT.0) THEN
         IERR = 2
         NREC = 1
         WRITE (REC,FMT=99998) M, N
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
         IERR = 3
         NREC = 2
         IF (NOTA) THEN
            WRITE (REC,FMT=99997) 'TRANSA', TRANSA, 'LDA', 'M', 'LDA',
     *        LDA, 'M', M
         ELSE
            WRITE (REC,FMT=99997) 'TRANSA', TRANSA, 'LDA', 'N', 'LDA',
     *        LDA, 'N', N
         END IF
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
         IERR = 4
         NREC = 2
         IF (NOTB) THEN
            WRITE (REC,FMT=99997) 'TRANSB', TRANSB, 'LDB', 'M', 'LDB',
     *        LDB, 'M', M
         ELSE
            WRITE (REC,FMT=99997) 'TRANSB', TRANSB, 'LDB', 'N', 'LDB',
     *        LDB, 'N', N
         END IF
      ELSE IF (LDC.LT.MAX(1,M)) THEN
         IERR = 5
         NREC = 1
         WRITE (REC,FMT=99996) LDC, M
      END IF
C
      IF (IERR.NE.0) THEN
         IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
         RETURN
      ELSE
         IFAIL = 0
      END IF
C
C     Quick return if possible.
C
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
C
C     Start the operations.
C
      IF (NOTA) THEN
         IF (NOTB) THEN
C
C           Form  C := alpha*A + beta*B.
C
            IF (ALPHA.EQ.ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 40 J = 1, N
                     DO 20 I = 1, M
                        C(I,J) = A(I,J) + B(I,J)
   20                CONTINUE
   40             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 80 J = 1, N
                     DO 60 I = 1, M
                        C(I,J) = A(I,J) - B(I,J)
   60                CONTINUE
   80             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 120 J = 1, N
                     DO 100 I = 1, M
                        C(I,J) = A(I,J)
  100                CONTINUE
  120             CONTINUE
               ELSE
                  DO 160 J = 1, N
                     DO 140 I = 1, M
                        C(I,J) = A(I,J) + BETA*B(I,J)
  140                CONTINUE
  160             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.-ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 200 J = 1, N
                     DO 180 I = 1, M
                        C(I,J) = -A(I,J) + B(I,J)
  180                CONTINUE
  200             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 240 J = 1, N
                     DO 220 I = 1, M
                        C(I,J) = -A(I,J) - B(I,J)
  220                CONTINUE
  240             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 280 J = 1, N
                     DO 260 I = 1, M
                        C(I,J) = -A(I,J)
  260                CONTINUE
  280             CONTINUE
               ELSE
                  DO 320 J = 1, N
                     DO 300 I = 1, M
                        C(I,J) = -A(I,J) + BETA*B(I,J)
  300                CONTINUE
  320             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.ZERO) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 360 J = 1, N
                     DO 340 I = 1, M
                        C(I,J) = B(I,J)
  340                CONTINUE
  360             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 400 J = 1, N
                     DO 380 I = 1, M
                        C(I,J) = -B(I,J)
  380                CONTINUE
  400             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 440 J = 1, N
                     DO 420 I = 1, M
                        C(I,J) = ZERO
  420                CONTINUE
  440             CONTINUE
               ELSE
                  DO 480 J = 1, N
                     DO 460 I = 1, M
                        C(I,J) = BETA*B(I,J)
  460                CONTINUE
  480             CONTINUE
               END IF
            ELSE
               IF (BETA.EQ.ONE) THEN
                  DO 520 J = 1, N
                     DO 500 I = 1, M
                        C(I,J) = ALPHA*A(I,J) + B(I,J)
  500                CONTINUE
  520             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 560 J = 1, N
                     DO 540 I = 1, M
                        C(I,J) = ALPHA*A(I,J) - B(I,J)
  540                CONTINUE
  560             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 600 J = 1, N
                     DO 580 I = 1, M
                        C(I,J) = ALPHA*A(I,J)
  580                CONTINUE
  600             CONTINUE
               ELSE
                  DO 640 J = 1, N
                     DO 620 I = 1, M
                        C(I,J) = ALPHA*A(I,J) + BETA*B(I,J)
  620                CONTINUE
  640             CONTINUE
               END IF
            END IF
         ELSE
C
C           Form  C := alpha*A + beta*B'.
C
            IF (ALPHA.EQ.ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 680 J = 1, N
                     DO 660 I = 1, M
                        C(I,J) = A(I,J) + B(J,I)
  660                CONTINUE
  680             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 720 J = 1, N
                     DO 700 I = 1, M
                        C(I,J) = A(I,J) - B(J,I)
  700                CONTINUE
  720             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 760 J = 1, N
                     DO 740 I = 1, M
                        C(I,J) = A(I,J)
  740                CONTINUE
  760             CONTINUE
               ELSE
                  DO 800 J = 1, N
                     DO 780 I = 1, M
                        C(I,J) = A(I,J) + BETA*B(J,I)
  780                CONTINUE
  800             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.-ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 840 J = 1, N
                     DO 820 I = 1, M
                        C(I,J) = -A(I,J) + B(J,I)
  820                CONTINUE
  840             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 880 J = 1, N
                     DO 860 I = 1, M
                        C(I,J) = -A(I,J) - B(J,I)
  860                CONTINUE
  880             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 920 J = 1, N
                     DO 900 I = 1, M
                        C(I,J) = -A(I,J)
  900                CONTINUE
  920             CONTINUE
               ELSE
                  DO 960 J = 1, N
                     DO 940 I = 1, M
                        C(I,J) = -A(I,J) + BETA*B(J,I)
  940                CONTINUE
  960             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.ZERO) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 1000 J = 1, N
                     DO 980 I = 1, M
                        C(I,J) = B(J,I)
  980                CONTINUE
 1000             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 1040 J = 1, N
                     DO 1020 I = 1, M
                        C(I,J) = -B(J,I)
 1020                CONTINUE
 1040             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 1080 J = 1, N
                     DO 1060 I = 1, M
                        C(I,J) = ZERO
 1060                CONTINUE
 1080             CONTINUE
               ELSE
                  DO 1120 J = 1, N
                     DO 1100 I = 1, M
                        C(I,J) = BETA*B(J,I)
 1100                CONTINUE
 1120             CONTINUE
               END IF
            ELSE
               IF (BETA.EQ.ONE) THEN
                  DO 1160 J = 1, N
                     DO 1140 I = 1, M
                        C(I,J) = ALPHA*A(I,J) + B(J,I)
 1140                CONTINUE
 1160             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 1200 J = 1, N
                     DO 1180 I = 1, M
                        C(I,J) = ALPHA*A(I,J) - B(J,I)
 1180                CONTINUE
 1200             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 1240 J = 1, N
                     DO 1220 I = 1, M
                        C(I,J) = ALPHA*A(I,J)
 1220                CONTINUE
 1240             CONTINUE
               ELSE
                  DO 1280 J = 1, N
                     DO 1260 I = 1, M
                        C(I,J) = ALPHA*A(I,J) + BETA*B(J,I)
 1260                CONTINUE
 1280             CONTINUE
               END IF
            END IF
         END IF
      ELSE
         IF (NOTB) THEN
C
C           Form  C := alpha*A' + beta*B.
C
            IF (ALPHA.EQ.ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 1320 J = 1, N
                     DO 1300 I = 1, M
                        C(I,J) = A(J,I) + B(I,J)
 1300                CONTINUE
 1320             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 1360 J = 1, N
                     DO 1340 I = 1, M
                        C(I,J) = A(J,I) - B(I,J)
 1340                CONTINUE
 1360             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 1400 J = 1, N
                     DO 1380 I = 1, M
                        C(I,J) = A(J,I)
 1380                CONTINUE
 1400             CONTINUE
               ELSE
                  DO 1440 J = 1, N
                     DO 1420 I = 1, M
                        C(I,J) = A(J,I) + BETA*B(I,J)
 1420                CONTINUE
 1440             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.-ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 1480 J = 1, N
                     DO 1460 I = 1, M
                        C(I,J) = -A(J,I) + B(I,J)
 1460                CONTINUE
 1480             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 1520 J = 1, N
                     DO 1500 I = 1, M
                        C(I,J) = -A(J,I) - B(I,J)
 1500                CONTINUE
 1520             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 1560 J = 1, N
                     DO 1540 I = 1, M
                        C(I,J) = -A(J,I)
 1540                CONTINUE
 1560             CONTINUE
               ELSE
                  DO 1600 J = 1, N
                     DO 1580 I = 1, M
                        C(I,J) = -A(J,I) + BETA*B(I,J)
 1580                CONTINUE
 1600             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.ZERO) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 1640 J = 1, N
                     DO 1620 I = 1, M
                        C(I,J) = B(I,J)
 1620                CONTINUE
 1640             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 1680 J = 1, N
                     DO 1660 I = 1, M
                        C(I,J) = -B(I,J)
 1660                CONTINUE
 1680             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 1720 J = 1, N
                     DO 1700 I = 1, M
                        C(I,J) = ZERO
 1700                CONTINUE
 1720             CONTINUE
               ELSE
                  DO 1760 J = 1, N
                     DO 1740 I = 1, M
                        C(I,J) = BETA*B(I,J)
 1740                CONTINUE
 1760             CONTINUE
               END IF
            ELSE
               IF (BETA.EQ.ONE) THEN
                  DO 1800 J = 1, N
                     DO 1780 I = 1, M
                        C(I,J) = ALPHA*A(J,I) + B(I,J)
 1780                CONTINUE
 1800             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 1840 J = 1, N
                     DO 1820 I = 1, M
                        C(I,J) = ALPHA*A(J,I) - B(I,J)
 1820                CONTINUE
 1840             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 1880 J = 1, N
                     DO 1860 I = 1, M
                        C(I,J) = ALPHA*A(J,I)
 1860                CONTINUE
 1880             CONTINUE
               ELSE
                  DO 1920 J = 1, N
                     DO 1900 I = 1, M
                        C(I,J) = ALPHA*A(J,I) + BETA*B(I,J)
 1900                CONTINUE
 1920             CONTINUE
               END IF
            END IF
         ELSE
C
C           Form  C := alpha*A' + beta*B.'
C
            IF (ALPHA.EQ.ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 1960 I = 1, M
                     DO 1940 J = 1, N
                        C(I,J) = A(J,I) + B(J,I)
 1940                CONTINUE
 1960             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 2000 I = 1, M
                     DO 1980 J = 1, N
                        C(I,J) = A(J,I) - B(J,I)
 1980                CONTINUE
 2000             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 2040 I = 1, M
                     DO 2020 J = 1, N
                        C(I,J) = A(J,I)
 2020                CONTINUE
 2040             CONTINUE
               ELSE
                  DO 2080 I = 1, M
                     DO 2060 J = 1, N
                        C(I,J) = A(J,I) + BETA*B(J,I)
 2060                CONTINUE
 2080             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.-ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 2120 I = 1, M
                     DO 2100 J = 1, N
                        C(I,J) = -A(J,I) + B(J,I)
 2100                CONTINUE
 2120             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 2160 I = 1, M
                     DO 2140 J = 1, N
                        C(I,J) = -A(J,I) - B(J,I)
 2140                CONTINUE
 2160             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 2200 I = 1, M
                     DO 2180 J = 1, N
                        C(I,J) = -A(J,I)
 2180                CONTINUE
 2200             CONTINUE
               ELSE
                  DO 2240 I = 1, M
                     DO 2220 J = 1, N
                        C(I,J) = -A(J,I) + BETA*B(J,I)
 2220                CONTINUE
 2240             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.ZERO) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 2280 I = 1, M
                     DO 2260 J = 1, N
                        C(I,J) = B(J,I)
 2260                CONTINUE
 2280             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 2320 I = 1, M
                     DO 2300 J = 1, N
                        C(I,J) = -B(J,I)
 2300                CONTINUE
 2320             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 2360 J = 1, N
                     DO 2340 I = 1, M
                        C(I,J) = ZERO
 2340                CONTINUE
 2360             CONTINUE
               ELSE
                  DO 2400 I = 1, M
                     DO 2380 J = 1, N
                        C(I,J) = BETA*B(J,I)
 2380                CONTINUE
 2400             CONTINUE
               END IF
            ELSE
               IF (BETA.EQ.ONE) THEN
                  DO 2440 I = 1, M
                     DO 2420 J = 1, N
                        C(I,J) = ALPHA*A(J,I) + B(J,I)
 2420                CONTINUE
 2440             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 2480 I = 1, M
                     DO 2460 J = 1, N
                        C(I,J) = ALPHA*A(J,I) - B(J,I)
 2460                CONTINUE
 2480             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 2520 I = 1, M
                     DO 2500 J = 1, N
                        C(I,J) = ALPHA*A(J,I)
 2500                CONTINUE
 2520             CONTINUE
               ELSE
                  DO 2560 I = 1, M
                     DO 2540 J = 1, N
                        C(I,J) = ALPHA*A(J,I) + BETA*B(J,I)
 2540                CONTINUE
 2560             CONTINUE
               END IF
            END IF
         END IF
      END IF
C
      RETURN
C
C
C     End of F01CTF .
C
99999 FORMAT (1X,'** On entry, ',A,' .ne. ''N'', ''T'' or ''C'': ',A,
     *       ' = ''',A,'''')
99998 FORMAT (1X,'** On entry, either M .lt. 0 or N .lt. 0: M =',I12,
     *       ', N =',I12)
99997 FORMAT (1X,'** On entry with ',A,' = ''',A,''', ',A,' .lt. max (',
     *       '1,',A,') :',/4X,A,' =',I10,', ',A,' =',I10)
99996 FORMAT (1X,'** On entry, LDC .lt. max (1, M) : LDC = ',I10,', M ',
     *       '= ',I10)
      END
