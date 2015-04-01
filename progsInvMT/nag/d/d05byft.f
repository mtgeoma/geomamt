      SUBROUTINE D05BYF(IORDER,IQ,LENFW,WT,SW,LDSW,WORK,LWK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     <<<<<<<<<<<<<<<<<<<<<<<<<     >>>>>>>>>>>>>>>>>>>>>>>>>>
C
C              A routine for computing Fractional weights
C                 of BDF formulae of orders 4, 5 and 6.
C
C     M.S. Derakhshan.
C     --------------------------------------------------------------
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This routine computes the square root convolution weights
C     of BDF reducible rules of orders 4 to 6. It also computes
C     the values of the starting fractional weights associated with
C     these BDF formulae.
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     --------------------------------------------------------------
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D05BYF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IORDER, IQ, LDSW, LENFW, LWK
C     .. Array Arguments ..
      DOUBLE PRECISION  SW(LDSW,0:2*IORDER-2), WORK(LWK), WT(0:LENFW-1)
C     .. Local Scalars ..
      INTEGER           INFO, IS, L1, L2, LIQ, LL, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(3)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D05BYH, D05BYZ
C     .. Executable Statements ..
      INFO = 0
      IS = 2*IORDER - 1
      LIQ = 2**(IQ+1)
      LL = 4*LIQ
C
      IF (IORDER.LT.4 .OR. IORDER.GT.6) THEN
         INFO = 1
         NREC = 1
         WRITE (P01REC,FMT=99999) IORDER
         GO TO 20
C
      ELSE IF (IQ.LT.0) THEN
         INFO = 1
         NREC = 1
         WRITE (P01REC,FMT=99998) IQ
         GO TO 20
C
      ELSE IF (LDSW.LT.(LIQ+IS)) THEN
         INFO = 1
         NREC = 3
         WRITE (P01REC,FMT=99997) LDSW, IQ, IORDER, LIQ + IS
         GO TO 20
C
      ELSE IF (LENFW.LT.(2*LIQ)) THEN
         INFO = 1
         NREC = 3
         WRITE (P01REC,FMT=99996) LENFW, IQ, 2*LIQ
         GO TO 20
C
      ELSE IF (LWK.LT.LL) THEN
         INFO = 1
         NREC = 3
         WRITE (P01REC,FMT=99995) LWK, IQ, LL
         GO TO 20
C
      END IF
C
      L1 = 1
      L2 = L1 + 2*LIQ
C
      CALL D05BYH(IQ,IORDER,LIQ,WT,WORK(L1),WORK(L2))
C
      CALL D05BYZ(WT(0),IS,IQ,LIQ,SW,LDSW,WORK(L1),WORK(L2))
C
C
   20 IFAIL = P01ABF(IFAIL,INFO,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, IORDER lt 4 or IORDER gt 6: IORDER = ',I16)
99998 FORMAT (' ** On entry, IQ lt 0: IQ = ',I16)
99997 FORMAT (' ** On entry, LDSW is too small: LDSW = ',I16,/' ** IQ ',
     *       '= ',I16,' and IORDER = ',I16,/' ** LDSW should be greate',
     *       'r than or equal to',I16)
99996 FORMAT (' ** On entry, LENFW is too small: LENFW =',I16,/' ** IQ',
     *       ' = ',I16,/' ** LENFW should be greater than or equal to',
     *       I16)
99995 FORMAT (' ** On entry, LWK is too small: LWK =',I16,/' ** IQ = ',
     *       I16,/' ** LWK should be greater than or equal to',I16)
      END
