      SUBROUTINE D02NMW(METH,ELCO,TESCO)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     OLD NAME DFODE
C
C-----------------------------------------------------------------------
C D02NMW IS CALLED BY THE INTEGRATOR ROUTINE TO SET COEFFICIENTS
C NEEDED THERE.  THE COEFFICIENTS FOR THE CURRENT METHOD, AS
C GIVEN BY THE VALUE OF METH, ARE SET FOR ALL ORDERS AND SAVED.
C THE MAXIMUM ORDER ASSUMED HERE IS 12 IF METH = 1 AND 5 IF METH = 2.
C (A SMALLER VALUE OF THE MAXIMUM ORDER IS ALSO ALLOWED.)
C D02NMW IS CALLED ONCE AT THE BEGINNING OF THE PROBLEM,
C AND IS NOT CALLED AGAIN UNLESS AND UNTIL METH IS CHANGED.
C
C THE ELCO ARRAY CONTAINS THE BASIC METHOD COEFFICIENTS.
C THE COEFFICIENTS EL(I), 1 .LE. I .LE. NQ+1, FOR THE METHOD OF
C ORDER NQ ARE STORED IN ELCO(I,NQ).  THEY ARE GIVEN BY A GENETRATING
C POLYNOMIAL, I.E.,
C     L(X) = EL(1) + EL(2)*X + ... + EL(NQ+1)*X**NQ.
C FOR THE IMPLICIT ADAMS METHODS, L(X) IS GIVEN BY
C     DL/DX = (X+1)*(X+2)*...*(X+NQ-1)/FACTORIAL(NQ-1),    L(-1) = 0.
C FOR THE BDF METHODS, L(X) IS GIVEN BY
C     L(X) = (X+1)*(X+2)* ... *(X+NQ)/K,
C WHERE         K = FACTORIAL(NQ)*(1 + 1/2 + ... + 1/NQ).
C
C THE TESCO ARRAY CONTAINS TEST CONSTANTS USED FOR THE
C LOCAL ERROR TEST AND THE SELECTION OF STEP SIZE AND/OR ORDER.
C AT ORDER NQ, TESCO(K,NQ) IS USED FOR THE SELECTION OF STEP
C SIZE AT ORDER NQ - 1 IF K = 1, AT ORDER NQ IF K = 2, AND AT ORDER
C NQ + 1 IF K = 3.
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           METH
C     .. Array Arguments ..
      DOUBLE PRECISION  ELCO(13,12), TESCO(3,12)
C     .. Local Scalars ..
      DOUBLE PRECISION  AGAMQ, FNQ, FNQM1, PINT, RAGQ, RQ1FAC, RQFAC,
     *                  TSIGN, XPIN
      INTEGER           I, IB, NQ, NQM1, NQP1
C     .. Local Arrays ..
      DOUBLE PRECISION  PC(12)
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
      GO TO (20,120) METH
C
   20 ELCO(1,1) = 1.0D0
      ELCO(2,1) = 1.0D0
      TESCO(1,1) = 0.0D0
      TESCO(2,1) = 2.0D0
      TESCO(1,2) = 1.0D0
      TESCO(3,12) = 0.0D0
      PC(1) = 1.0D0
      RQFAC = 1.0D0
      DO 100 NQ = 2, 12
C-----------------------------------------------------------------------
C THE PC ARRAY WILL CONTAIN THE COEFFICIENTS OF THE POLYNOMIAL
C     P(X) = (X+1)*(X+2)*...*(X+NQ-1).
C INITIALLY, P(X) = 1.
C-----------------------------------------------------------------------
         RQ1FAC = RQFAC
         RQFAC = RQFAC/DBLE(NQ)
         NQM1 = NQ - 1
         FNQM1 = DBLE(NQM1)
         NQP1 = NQ + 1
C FORM COEFFICIENTS OF P(X)*(X+NQ-1). ----------------------------------
         PC(NQ) = 0.0D0
         DO 40 IB = 1, NQM1
            I = NQP1 - IB
            PC(I) = PC(I-1) + FNQM1*PC(I)
   40    CONTINUE
         PC(1) = FNQM1*PC(1)
C COMPUTE INTEGRAL, -1 TO 0, OF P(X) AND X*P(X). -----------------------
         PINT = PC(1)
         XPIN = PC(1)/2.0D0
         TSIGN = 1.0D0
         DO 60 I = 2, NQ
            TSIGN = -TSIGN
            PINT = PINT + TSIGN*PC(I)/DBLE(I)
            XPIN = XPIN + TSIGN*PC(I)/DBLE(I+1)
   60    CONTINUE
C STORE COEFFICIENTS IN ELCO AND TESCO. --------------------------------
         ELCO(1,NQ) = PINT*RQ1FAC
         ELCO(2,NQ) = 1.0D0
         DO 80 I = 2, NQ
            ELCO(I+1,NQ) = RQ1FAC*PC(I)/DBLE(I)
   80    CONTINUE
         AGAMQ = RQFAC*XPIN
         RAGQ = 1.0D0/AGAMQ
         TESCO(2,NQ) = RAGQ
         IF (NQ.LT.12) TESCO(1,NQP1) = RAGQ*RQFAC/DBLE(NQP1)
         TESCO(3,NQM1) = RAGQ
  100 CONTINUE
      RETURN
C
  120 PC(1) = 1.0D0
      RQ1FAC = 1.0D0
      DO 180 NQ = 1, 5
C-----------------------------------------------------------------------
C THE PC ARRAY WILL CONTAIN THE COEFFICIENTS OF THE POLYNOMIAL
C     P(X) = (X+1)*(X+2)*...*(X+NQ).
C INITIALLY, P(X) = 1.
C-----------------------------------------------------------------------
         FNQ = DBLE(NQ)
         NQP1 = NQ + 1
C FORM COEFFICIENTS OF P(X)*(X+NQ). ------------------------------------
         PC(NQP1) = 0.0D0
         DO 140 IB = 1, NQ
            I = NQ + 2 - IB
            PC(I) = PC(I-1) + FNQ*PC(I)
  140    CONTINUE
         PC(1) = FNQ*PC(1)
C STORE COEFFICIENTS IN ELCO AND TESCO. --------------------------------
         DO 160 I = 1, NQP1
            ELCO(I,NQ) = PC(I)/PC(2)
  160    CONTINUE
         ELCO(2,NQ) = 1.0D0
         TESCO(1,NQ) = RQ1FAC
         TESCO(2,NQ) = DBLE(NQP1)/ELCO(1,NQ)
         TESCO(3,NQ) = DBLE(NQ+2)/ELCO(1,NQ)
         RQ1FAC = RQ1FAC/FNQ
  180 CONTINUE
      RETURN
      END
