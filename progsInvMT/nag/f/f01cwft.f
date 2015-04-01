      SUBROUTINE F01CWF(TRANSA,TRANSB,M,N,ALPHA,A,LDA,BETA,B,LDB,C,LDC,
     *                  IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15A REVISED. IER-909 (APR 1991).
C
C  Purpose
C  =======
C
C  F01CWF performs one of the matrix-matrix operations
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
C              TRANSA = 'C' or 'c',  op( A ) = conjg(A').
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
C              TRANSB = 'C' or 'c',  op( B ) = conjg(B').
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
C  ALPHA  - COMPLEX         .
C           On entry, ALPHA specifies the scalar alpha.  When  ALPHA  is
C           supplied as zero then A need not be set on entry.
C           Unchanged on exit.
C
C  A      - COMPLEX          array of DIMENSION ( LDA, ma ), where ma is
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
C  BETA   - COMPLEX         .
C           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
C           supplied as zero then B need not be set on entry.
C           Unchanged on exit.
C
C  B      - COMPLEX          array of DIMENSION ( LDB, mb ), where mb is
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
C  C      - COMPLEX          array of DIMENSION ( LDC, N ).
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
      COMPLEX*16        ONE, ZERO
      PARAMETER         (ONE=(1.0D+0,0.0D+0),ZERO=(0.0D+0,0.0D+0))
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01CWF')
C     .. Scalar Arguments ..
      COMPLEX*16        ALPHA, BETA
      INTEGER           IFAIL, LDA, LDB, LDC, M, N
      CHARACTER*1       TRANSA, TRANSB
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*), C(LDC,N)
C     .. Local Scalars ..
      INTEGER           I, IERR, J, NREC, NROWA, NROWB
      LOGICAL           CONJA, CONJB, NOTA, NOTB
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG, MAX
C     .. Executable Statements ..
C
C     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
C     transposed.
C
      NOTA = TRANSA .EQ. 'N' .OR. TRANSA .EQ. 'n'
      NOTB = TRANSB .EQ. 'N' .OR. TRANSB .EQ. 'n'
      CONJA = TRANSA .EQ. 'C' .OR. TRANSA .EQ. 'c'
      CONJB = TRANSB .EQ. 'C' .OR. TRANSB .EQ. 'c'
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
      IF ( .NOT. NOTA .AND. .NOT. CONJA .AND. .NOT.
     *    (TRANSA.EQ.'T' .OR. TRANSA.EQ.'t')) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99999) 'TRANSA', 'TRANSA', TRANSA
      ELSE IF ( .NOT. NOTB .AND. .NOT. CONJB .AND. .NOT.
     *         (TRANSB.EQ.'T' .OR. TRANSB.EQ.'t')) THEN
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
         ELSE IF (CONJB) THEN
C
C           Form  C := alpha*A + beta*conjg(B').
C
            IF (ALPHA.EQ.ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 680 J = 1, N
                     DO 660 I = 1, M
                        C(I,J) = A(I,J) + DCONJG(B(J,I))
  660                CONTINUE
  680             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 720 J = 1, N
                     DO 700 I = 1, M
                        C(I,J) = A(I,J) - DCONJG(B(J,I))
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
                        C(I,J) = A(I,J) + BETA*DCONJG(B(J,I))
  780                CONTINUE
  800             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.-ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 840 J = 1, N
                     DO 820 I = 1, M
                        C(I,J) = -A(I,J) + DCONJG(B(J,I))
  820                CONTINUE
  840             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 880 J = 1, N
                     DO 860 I = 1, M
                        C(I,J) = -A(I,J) - DCONJG(B(J,I))
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
                        C(I,J) = -A(I,J) + BETA*DCONJG(B(J,I))
  940                CONTINUE
  960             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.ZERO) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 1000 J = 1, N
                     DO 980 I = 1, M
                        C(I,J) = DCONJG(B(J,I))
  980                CONTINUE
 1000             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 1040 J = 1, N
                     DO 1020 I = 1, M
                        C(I,J) = -DCONJG(B(J,I))
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
                        C(I,J) = BETA*DCONJG(B(J,I))
 1100                CONTINUE
 1120             CONTINUE
               END IF
            ELSE
               IF (BETA.EQ.ONE) THEN
                  DO 1160 J = 1, N
                     DO 1140 I = 1, M
                        C(I,J) = ALPHA*A(I,J) + DCONJG(B(J,I))
 1140                CONTINUE
 1160             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 1200 J = 1, N
                     DO 1180 I = 1, M
                        C(I,J) = ALPHA*A(I,J) - DCONJG(B(J,I))
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
                        C(I,J) = ALPHA*A(I,J) + BETA*DCONJG(B(J,I))
 1260                CONTINUE
 1280             CONTINUE
               END IF
            END IF
         ELSE
C
C           Form  C := alpha*A + beta*B'.
C
            IF (ALPHA.EQ.ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 1320 J = 1, N
                     DO 1300 I = 1, M
                        C(I,J) = A(I,J) + B(J,I)
 1300                CONTINUE
 1320             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 1360 J = 1, N
                     DO 1340 I = 1, M
                        C(I,J) = A(I,J) - B(J,I)
 1340                CONTINUE
 1360             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 1400 J = 1, N
                     DO 1380 I = 1, M
                        C(I,J) = A(I,J)
 1380                CONTINUE
 1400             CONTINUE
               ELSE
                  DO 1440 J = 1, N
                     DO 1420 I = 1, M
                        C(I,J) = A(I,J) + BETA*B(J,I)
 1420                CONTINUE
 1440             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.-ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 1480 J = 1, N
                     DO 1460 I = 1, M
                        C(I,J) = -A(I,J) + B(J,I)
 1460                CONTINUE
 1480             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 1520 J = 1, N
                     DO 1500 I = 1, M
                        C(I,J) = -A(I,J) - B(J,I)
 1500                CONTINUE
 1520             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 1560 J = 1, N
                     DO 1540 I = 1, M
                        C(I,J) = -A(I,J)
 1540                CONTINUE
 1560             CONTINUE
               ELSE
                  DO 1600 J = 1, N
                     DO 1580 I = 1, M
                        C(I,J) = -A(I,J) + BETA*B(J,I)
 1580                CONTINUE
 1600             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.ZERO) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 1640 J = 1, N
                     DO 1620 I = 1, M
                        C(I,J) = B(J,I)
 1620                CONTINUE
 1640             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 1680 J = 1, N
                     DO 1660 I = 1, M
                        C(I,J) = -B(J,I)
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
                        C(I,J) = BETA*B(J,I)
 1740                CONTINUE
 1760             CONTINUE
               END IF
            ELSE
               IF (BETA.EQ.ONE) THEN
                  DO 1800 J = 1, N
                     DO 1780 I = 1, M
                        C(I,J) = ALPHA*A(I,J) + B(J,I)
 1780                CONTINUE
 1800             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 1840 J = 1, N
                     DO 1820 I = 1, M
                        C(I,J) = ALPHA*A(I,J) - B(J,I)
 1820                CONTINUE
 1840             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 1880 J = 1, N
                     DO 1860 I = 1, M
                        C(I,J) = ALPHA*A(I,J)
 1860                CONTINUE
 1880             CONTINUE
               ELSE
                  DO 1920 J = 1, N
                     DO 1900 I = 1, M
                        C(I,J) = ALPHA*A(I,J) + BETA*B(J,I)
 1900                CONTINUE
 1920             CONTINUE
               END IF
            END IF
         END IF
      ELSE IF (CONJA) THEN
         IF (NOTB) THEN
C
C           Form  C := alpha*conjg(A') + beta*B.
C
            IF (ALPHA.EQ.ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 1960 J = 1, N
                     DO 1940 I = 1, M
                        C(I,J) = DCONJG(A(J,I)) + B(I,J)
 1940                CONTINUE
 1960             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 2000 J = 1, N
                     DO 1980 I = 1, M
                        C(I,J) = DCONJG(A(J,I)) - B(I,J)
 1980                CONTINUE
 2000             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 2040 J = 1, N
                     DO 2020 I = 1, M
                        C(I,J) = DCONJG(A(J,I))
 2020                CONTINUE
 2040             CONTINUE
               ELSE
                  DO 2080 J = 1, N
                     DO 2060 I = 1, M
                        C(I,J) = DCONJG(A(J,I)) + BETA*B(I,J)
 2060                CONTINUE
 2080             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.-ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 2120 J = 1, N
                     DO 2100 I = 1, M
                        C(I,J) = -DCONJG(A(J,I)) + B(I,J)
 2100                CONTINUE
 2120             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 2160 J = 1, N
                     DO 2140 I = 1, M
                        C(I,J) = -DCONJG(A(J,I)) - B(I,J)
 2140                CONTINUE
 2160             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 2200 J = 1, N
                     DO 2180 I = 1, M
                        C(I,J) = -DCONJG(A(J,I))
 2180                CONTINUE
 2200             CONTINUE
               ELSE
                  DO 2240 J = 1, N
                     DO 2220 I = 1, M
                        C(I,J) = -DCONJG(A(J,I)) + BETA*B(I,J)
 2220                CONTINUE
 2240             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.ZERO) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 2280 J = 1, N
                     DO 2260 I = 1, M
                        C(I,J) = B(I,J)
 2260                CONTINUE
 2280             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 2320 J = 1, N
                     DO 2300 I = 1, M
                        C(I,J) = -B(I,J)
 2300                CONTINUE
 2320             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 2360 J = 1, N
                     DO 2340 I = 1, M
                        C(I,J) = ZERO
 2340                CONTINUE
 2360             CONTINUE
               ELSE
                  DO 2400 J = 1, N
                     DO 2380 I = 1, M
                        C(I,J) = BETA*B(I,J)
 2380                CONTINUE
 2400             CONTINUE
               END IF
            ELSE
               IF (BETA.EQ.ONE) THEN
                  DO 2440 J = 1, N
                     DO 2420 I = 1, M
                        C(I,J) = ALPHA*DCONJG(A(J,I)) + B(I,J)
 2420                CONTINUE
 2440             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 2480 J = 1, N
                     DO 2460 I = 1, M
                        C(I,J) = ALPHA*DCONJG(A(J,I)) - B(I,J)
 2460                CONTINUE
 2480             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 2520 J = 1, N
                     DO 2500 I = 1, M
                        C(I,J) = ALPHA*DCONJG(A(J,I))
 2500                CONTINUE
 2520             CONTINUE
               ELSE
                  DO 2560 J = 1, N
                     DO 2540 I = 1, M
                        C(I,J) = ALPHA*DCONJG(A(J,I)) + BETA*B(I,J)
 2540                CONTINUE
 2560             CONTINUE
               END IF
            END IF
         ELSE IF (CONJB) THEN
C
C           Form  C := alpha*conjg(A') + beta*conjg(B').
C
            IF (ALPHA.EQ.ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 2600 I = 1, M
                     DO 2580 J = 1, N
                        C(I,J) = DCONJG(A(J,I)) + DCONJG(B(J,I))
 2580                CONTINUE
 2600             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 2640 I = 1, M
                     DO 2620 J = 1, N
                        C(I,J) = DCONJG(A(J,I)) - DCONJG(B(J,I))
 2620                CONTINUE
 2640             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 2680 I = 1, M
                     DO 2660 J = 1, N
                        C(I,J) = DCONJG(A(J,I))
 2660                CONTINUE
 2680             CONTINUE
               ELSE
                  DO 2720 I = 1, M
                     DO 2700 J = 1, N
                        C(I,J) = DCONJG(A(J,I)) + BETA*DCONJG(B(J,I))
 2700                CONTINUE
 2720             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.-ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 2760 I = 1, M
                     DO 2740 J = 1, N
                        C(I,J) = -DCONJG(A(J,I)) + DCONJG(B(J,I))
 2740                CONTINUE
 2760             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 2800 I = 1, M
                     DO 2780 J = 1, N
                        C(I,J) = -DCONJG(A(J,I)) - DCONJG(B(J,I))
 2780                CONTINUE
 2800             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 2840 I = 1, M
                     DO 2820 J = 1, N
                        C(I,J) = -DCONJG(A(J,I))
 2820                CONTINUE
 2840             CONTINUE
               ELSE
                  DO 2880 I = 1, M
                     DO 2860 J = 1, N
                        C(I,J) = -DCONJG(A(J,I)) + BETA*DCONJG(B(J,I))
 2860                CONTINUE
 2880             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.ZERO) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 2920 I = 1, M
                     DO 2900 J = 1, N
                        C(I,J) = DCONJG(B(J,I))
 2900                CONTINUE
 2920             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 2960 I = 1, M
                     DO 2940 J = 1, N
                        C(I,J) = -DCONJG(B(J,I))
 2940                CONTINUE
 2960             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 3000 J = 1, N
                     DO 2980 I = 1, M
                        C(I,J) = ZERO
 2980                CONTINUE
 3000             CONTINUE
               ELSE
                  DO 3040 I = 1, M
                     DO 3020 J = 1, N
                        C(I,J) = BETA*DCONJG(B(J,I))
 3020                CONTINUE
 3040             CONTINUE
               END IF
            ELSE
               IF (BETA.EQ.ONE) THEN
                  DO 3080 I = 1, M
                     DO 3060 J = 1, N
                        C(I,J) = ALPHA*DCONJG(A(J,I)) + DCONJG(B(J,I))
 3060                CONTINUE
 3080             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 3120 I = 1, M
                     DO 3100 J = 1, N
                        C(I,J) = ALPHA*DCONJG(A(J,I)) - DCONJG(B(J,I))
 3100                CONTINUE
 3120             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 3160 I = 1, M
                     DO 3140 J = 1, N
                        C(I,J) = ALPHA*DCONJG(A(J,I))
 3140                CONTINUE
 3160             CONTINUE
               ELSE
                  DO 3200 I = 1, M
                     DO 3180 J = 1, N
                        C(I,J) = ALPHA*DCONJG(A(J,I)) +
     *                           BETA*DCONJG(B(J,I))
 3180                CONTINUE
 3200             CONTINUE
               END IF
            END IF
         ELSE
C
C           Form  C := alpha*conjg(A') + beta*B'.
C
            IF (ALPHA.EQ.ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 3240 I = 1, M
                     DO 3220 J = 1, N
                        C(I,J) = DCONJG(A(J,I)) + B(J,I)
 3220                CONTINUE
 3240             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 3280 I = 1, M
                     DO 3260 J = 1, N
                        C(I,J) = DCONJG(A(J,I)) - B(J,I)
 3260                CONTINUE
 3280             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 3320 I = 1, M
                     DO 3300 J = 1, N
                        C(I,J) = DCONJG(A(J,I))
 3300                CONTINUE
 3320             CONTINUE
               ELSE
                  DO 3360 I = 1, M
                     DO 3340 J = 1, N
                        C(I,J) = DCONJG(A(J,I)) + BETA*B(J,I)
 3340                CONTINUE
 3360             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.-ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 3400 I = 1, M
                     DO 3380 J = 1, N
                        C(I,J) = -DCONJG(A(J,I)) + B(J,I)
 3380                CONTINUE
 3400             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 3440 I = 1, M
                     DO 3420 J = 1, N
                        C(I,J) = -DCONJG(A(J,I)) - B(J,I)
 3420                CONTINUE
 3440             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 3480 I = 1, M
                     DO 3460 J = 1, N
                        C(I,J) = -DCONJG(A(J,I))
 3460                CONTINUE
 3480             CONTINUE
               ELSE
                  DO 3520 I = 1, M
                     DO 3500 J = 1, N
                        C(I,J) = -DCONJG(A(J,I)) + BETA*B(J,I)
 3500                CONTINUE
 3520             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.ZERO) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 3560 I = 1, M
                     DO 3540 J = 1, N
                        C(I,J) = B(J,I)
 3540                CONTINUE
 3560             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 3600 I = 1, M
                     DO 3580 J = 1, N
                        C(I,J) = -B(J,I)
 3580                CONTINUE
 3600             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 3640 J = 1, N
                     DO 3620 I = 1, M
                        C(I,J) = ZERO
 3620                CONTINUE
 3640             CONTINUE
               ELSE
                  DO 3680 I = 1, M
                     DO 3660 J = 1, N
                        C(I,J) = BETA*B(J,I)
 3660                CONTINUE
 3680             CONTINUE
               END IF
            ELSE
               IF (BETA.EQ.ONE) THEN
                  DO 3720 I = 1, M
                     DO 3700 J = 1, N
                        C(I,J) = ALPHA*DCONJG(A(J,I)) + B(J,I)
 3700                CONTINUE
 3720             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 3760 I = 1, M
                     DO 3740 J = 1, N
                        C(I,J) = ALPHA*DCONJG(A(J,I)) - B(J,I)
 3740                CONTINUE
 3760             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 3800 I = 1, M
                     DO 3780 J = 1, N
                        C(I,J) = ALPHA*DCONJG(A(J,I))
 3780                CONTINUE
 3800             CONTINUE
               ELSE
                  DO 3840 I = 1, M
                     DO 3820 J = 1, N
                        C(I,J) = ALPHA*DCONJG(A(J,I)) + BETA*B(J,I)
 3820                CONTINUE
 3840             CONTINUE
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
                  DO 3880 J = 1, N
                     DO 3860 I = 1, M
                        C(I,J) = A(J,I) + B(I,J)
 3860                CONTINUE
 3880             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 3920 J = 1, N
                     DO 3900 I = 1, M
                        C(I,J) = A(J,I) - B(I,J)
 3900                CONTINUE
 3920             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 3960 J = 1, N
                     DO 3940 I = 1, M
                        C(I,J) = A(J,I)
 3940                CONTINUE
 3960             CONTINUE
               ELSE
                  DO 4000 J = 1, N
                     DO 3980 I = 1, M
                        C(I,J) = A(J,I) + BETA*B(I,J)
 3980                CONTINUE
 4000             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.-ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 4040 J = 1, N
                     DO 4020 I = 1, M
                        C(I,J) = -A(J,I) + B(I,J)
 4020                CONTINUE
 4040             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 4080 J = 1, N
                     DO 4060 I = 1, M
                        C(I,J) = -A(J,I) - B(I,J)
 4060                CONTINUE
 4080             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 4120 J = 1, N
                     DO 4100 I = 1, M
                        C(I,J) = -A(J,I)
 4100                CONTINUE
 4120             CONTINUE
               ELSE
                  DO 4160 J = 1, N
                     DO 4140 I = 1, M
                        C(I,J) = -A(J,I) + BETA*B(I,J)
 4140                CONTINUE
 4160             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.ZERO) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 4200 J = 1, N
                     DO 4180 I = 1, M
                        C(I,J) = B(I,J)
 4180                CONTINUE
 4200             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 4240 J = 1, N
                     DO 4220 I = 1, M
                        C(I,J) = -B(I,J)
 4220                CONTINUE
 4240             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 4280 J = 1, N
                     DO 4260 I = 1, M
                        C(I,J) = ZERO
 4260                CONTINUE
 4280             CONTINUE
               ELSE
                  DO 4320 J = 1, N
                     DO 4300 I = 1, M
                        C(I,J) = BETA*B(I,J)
 4300                CONTINUE
 4320             CONTINUE
               END IF
            ELSE
               IF (BETA.EQ.ONE) THEN
                  DO 4360 J = 1, N
                     DO 4340 I = 1, M
                        C(I,J) = ALPHA*A(J,I) + B(I,J)
 4340                CONTINUE
 4360             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 4400 J = 1, N
                     DO 4380 I = 1, M
                        C(I,J) = ALPHA*A(J,I) - B(I,J)
 4380                CONTINUE
 4400             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 4440 J = 1, N
                     DO 4420 I = 1, M
                        C(I,J) = ALPHA*A(J,I)
 4420                CONTINUE
 4440             CONTINUE
               ELSE
                  DO 4480 J = 1, N
                     DO 4460 I = 1, M
                        C(I,J) = ALPHA*A(J,I) + BETA*B(I,J)
 4460                CONTINUE
 4480             CONTINUE
               END IF
            END IF
         ELSE IF (CONJB) THEN
C
C           Form  C := alpha*A' + beta*conjg(B').
C
            IF (ALPHA.EQ.ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 4520 I = 1, M
                     DO 4500 J = 1, N
                        C(I,J) = A(J,I) + DCONJG(B(J,I))
 4500                CONTINUE
 4520             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 4560 I = 1, M
                     DO 4540 J = 1, N
                        C(I,J) = A(J,I) - DCONJG(B(J,I))
 4540                CONTINUE
 4560             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 4600 I = 1, M
                     DO 4580 J = 1, N
                        C(I,J) = A(J,I)
 4580                CONTINUE
 4600             CONTINUE
               ELSE
                  DO 4640 I = 1, M
                     DO 4620 J = 1, N
                        C(I,J) = A(J,I) + BETA*DCONJG(B(J,I))
 4620                CONTINUE
 4640             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.-ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 4680 I = 1, M
                     DO 4660 J = 1, N
                        C(I,J) = -A(J,I) + DCONJG(B(J,I))
 4660                CONTINUE
 4680             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 4720 I = 1, M
                     DO 4700 J = 1, N
                        C(I,J) = -A(J,I) - DCONJG(B(J,I))
 4700                CONTINUE
 4720             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 4760 I = 1, M
                     DO 4740 J = 1, N
                        C(I,J) = -A(J,I)
 4740                CONTINUE
 4760             CONTINUE
               ELSE
                  DO 4800 I = 1, M
                     DO 4780 J = 1, N
                        C(I,J) = -A(J,I) + BETA*DCONJG(B(J,I))
 4780                CONTINUE
 4800             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.ZERO) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 4840 I = 1, M
                     DO 4820 J = 1, N
                        C(I,J) = DCONJG(B(J,I))
 4820                CONTINUE
 4840             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 4880 I = 1, M
                     DO 4860 J = 1, N
                        C(I,J) = -DCONJG(B(J,I))
 4860                CONTINUE
 4880             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 4920 J = 1, N
                     DO 4900 I = 1, M
                        C(I,J) = ZERO
 4900                CONTINUE
 4920             CONTINUE
               ELSE
                  DO 4960 I = 1, M
                     DO 4940 J = 1, N
                        C(I,J) = BETA*DCONJG(B(J,I))
 4940                CONTINUE
 4960             CONTINUE
               END IF
            ELSE
               IF (BETA.EQ.ONE) THEN
                  DO 5000 I = 1, M
                     DO 4980 J = 1, N
                        C(I,J) = ALPHA*A(J,I) + DCONJG(B(J,I))
 4980                CONTINUE
 5000             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 5040 I = 1, M
                     DO 5020 J = 1, N
                        C(I,J) = ALPHA*A(J,I) - DCONJG(B(J,I))
 5020                CONTINUE
 5040             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 5080 I = 1, M
                     DO 5060 J = 1, N
                        C(I,J) = ALPHA*A(J,I)
 5060                CONTINUE
 5080             CONTINUE
               ELSE
                  DO 5120 I = 1, M
                     DO 5100 J = 1, N
                        C(I,J) = ALPHA*A(J,I) + BETA*DCONJG(B(J,I))
 5100                CONTINUE
 5120             CONTINUE
               END IF
            END IF
         ELSE
C
C           Form  C := alpha*A' + beta*B'.
C
            IF (ALPHA.EQ.ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 5160 I = 1, M
                     DO 5140 J = 1, N
                        C(I,J) = A(J,I) + B(J,I)
 5140                CONTINUE
 5160             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 5200 I = 1, M
                     DO 5180 J = 1, N
                        C(I,J) = A(J,I) - B(J,I)
 5180                CONTINUE
 5200             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 5240 I = 1, M
                     DO 5220 J = 1, N
                        C(I,J) = A(J,I)
 5220                CONTINUE
 5240             CONTINUE
               ELSE
                  DO 5280 I = 1, M
                     DO 5260 J = 1, N
                        C(I,J) = A(J,I) + BETA*B(J,I)
 5260                CONTINUE
 5280             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.-ONE) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 5320 I = 1, M
                     DO 5300 J = 1, N
                        C(I,J) = -A(J,I) + B(J,I)
 5300                CONTINUE
 5320             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 5360 I = 1, M
                     DO 5340 J = 1, N
                        C(I,J) = -A(J,I) - B(J,I)
 5340                CONTINUE
 5360             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 5400 I = 1, M
                     DO 5380 J = 1, N
                        C(I,J) = -A(J,I)
 5380                CONTINUE
 5400             CONTINUE
               ELSE
                  DO 5440 I = 1, M
                     DO 5420 J = 1, N
                        C(I,J) = -A(J,I) + BETA*B(J,I)
 5420                CONTINUE
 5440             CONTINUE
               END IF
            ELSE IF (ALPHA.EQ.ZERO) THEN
               IF (BETA.EQ.ONE) THEN
                  DO 5480 I = 1, M
                     DO 5460 J = 1, N
                        C(I,J) = B(J,I)
 5460                CONTINUE
 5480             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 5520 I = 1, M
                     DO 5500 J = 1, N
                        C(I,J) = -B(J,I)
 5500                CONTINUE
 5520             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 5560 J = 1, N
                     DO 5540 I = 1, M
                        C(I,J) = ZERO
 5540                CONTINUE
 5560             CONTINUE
               ELSE
                  DO 5600 I = 1, M
                     DO 5580 J = 1, N
                        C(I,J) = BETA*B(J,I)
 5580                CONTINUE
 5600             CONTINUE
               END IF
            ELSE
               IF (BETA.EQ.ONE) THEN
                  DO 5640 I = 1, M
                     DO 5620 J = 1, N
                        C(I,J) = ALPHA*A(J,I) + B(J,I)
 5620                CONTINUE
 5640             CONTINUE
               ELSE IF (BETA.EQ.-ONE) THEN
                  DO 5680 I = 1, M
                     DO 5660 J = 1, N
                        C(I,J) = ALPHA*A(J,I) - B(J,I)
 5660                CONTINUE
 5680             CONTINUE
               ELSE IF (BETA.EQ.ZERO) THEN
                  DO 5720 I = 1, M
                     DO 5700 J = 1, N
                        C(I,J) = ALPHA*A(J,I)
 5700                CONTINUE
 5720             CONTINUE
               ELSE
                  DO 5760 I = 1, M
                     DO 5740 J = 1, N
                        C(I,J) = ALPHA*A(J,I) + BETA*B(J,I)
 5740                CONTINUE
 5760             CONTINUE
               END IF
            END IF
         END IF
      END IF
C
      RETURN
C
C
C     End of F01CWF .
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
