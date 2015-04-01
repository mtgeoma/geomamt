      SUBROUTINE E01DAS(N,IBANDW,UFCTR,LUFCTR,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     A DASL routine, unmodified except for name.
C
C     **********************************************************
C
C     D A S L  -  DATA APPROXIMATION SUBROUTINE LIBRARY
C
C     SUBROUTINE CHKRB     CHECK INVERTIBILITY OF BAND
C     ================     UPPER TRIANGULAR MATRIX.
C
C     CREATED 07 07 80.  UPDATED 24 06 82.  RELEASE 00/06
C
C     AUTHORS ... MAURICE G. COX AND PAULINE E. M. CURTIS.
C     NATIONAL PHYSICAL LABORATORY, TEDDINGTON,
C     MIDDLESEX TW11 OLW, ENGLAND.
C
C     (C)  CROWN COPYRIGHT 1980-1982
C
C     **********************************************************
C
C     CHKRB.  AN ALGORITHM FOR DETERMINING WHETHER A BAND
C     UPPER TRIANGULAR MATRIX IS INVERTIBLE
C     BY EXAMINING ITS DIAGONAL ELEMENTS.
C
C     INPUT PARAMETERS
C        N        ORDER OF BAND UPPER TRIANGULAR MATRIX
C        IBANDW   BANDWIDTH OF BAND UPPER TRIANGULAR MATRIX
C        UFCTR    (BAND PART ONLY OF) BAND UPPER TRIANGULAR
C                    MATRIX, STORED CONTIGUOUSLY IN ROW
C                    ORDER
C        LUFCTR   DIMENSION OF  UFCTR
C
C     FAILURE INDICATOR PARAMETER
C        IFAIL    FAILURE INDICATOR
C                    1 - ONE OF THE DIAGONAL ELEMENTS OF
C                        UFCTR  IS ZERO
C
C     ----------------------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER           IBANDW, IFAIL, LUFCTR, N
C     .. Array Arguments ..
      DOUBLE PRECISION  UFCTR(LUFCTR)
C     .. Local Scalars ..
      DOUBLE PRECISION  ZERO
      INTEGER           IERROR, IU, J, JTEST, NP2
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Data statements ..
C
      DATA              ZERO/0.0D+0/
C     .. Executable Statements ..
C
      IERROR = 1
      NP2 = N + 2
      JTEST = NP2 - IBANDW
      IU = 1 - MIN(IBANDW,N+1)
      DO 20 J = 1, N
         IF (J.LE.JTEST) IU = IU + IBANDW
         IF (J.GT.JTEST) IU = IU + NP2 - J
         IF (UFCTR(IU).EQ.ZERO) GO TO 40
   20 CONTINUE
      IERROR = 0
   40 IFAIL = IERROR
      RETURN
C
C     END E01DAS
C
      END
