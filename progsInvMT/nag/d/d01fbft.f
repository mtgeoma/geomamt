      DOUBLE PRECISION FUNCTION D01FBF(NDIM,NPTVEC,LWA,WEIGHT,ABSCIS,
     *                                 FUN,IFAIL)
C     ***   MULTIDIMENSIONAL INTEGRATION RULE EVALUATION
C     **************************************************************
C
C       INPUT PARAMETERS
C     NDIM  INTEGER WHICH SPECIFIES THE NUMBER OF VARIABLES
C           NDIM MUST BE GREATER THAN ONE AND LESS THAN 21.
C     NPTVEC INTEGER ARRAY OF DIMENSION NDIM. NPTVEC(J) SPECIFIES
C           THE NUMBER OF FUNCTION EVALUATION POINTS IN THE
C           INTEGRATION RULE FOR VARIABLE J.
C     LWA   INTEGER WHICH SPECIFIES THE DIMENSION OF THE ARRAYS
C           WEIGHT AND ABSCIS. LWA .GE. THE SUM OF THE COMPONENTS
C           OF NPTVEC.
C     WEIGHT  REAL ARRAY OF DIMENSION LWA.  THIS ARRAY MUST
C           CONTAIN ALL THE INTEGRATION RULE WEIGHTS IN
C           THE CORRECT ORDER.
C     ABSCIS REAL ARRAY OF DIMENSION LWA. THIS ARRAY MUST
C           CONTAIN ALL THE INTEGRATION RULE ABSCISSAE IN
C           THE CORRECT ORDER.
C     FUN    EXTERNALLY DECLARED USER DEFINED INTEGRAND
C           WHICH MUST HAVE PARAMETERS NDIM AND A REAL
C           ARRAY OF DIMENSION NDIM.
C     IFAIL INTEGER NAG FAILURE PARAMETER.  SEE NAG DOCUMENTATION.
C           ON ENTRY IFAIL SHOULD ALWAYS BE SET TO ZERO
C
C       OUTPUT PARAMETERS
C     D01FBF REAL VALUE WHICH GIVES THE VALUE OF THE
C           INTEGRATION RULE WHEN APPLIED TO THE INTEGRAND
C           FUN.
C     IFAIL  INTEGER NAG FAILURE PARAMETER
C           IFAIL=0 FOR NORMAL EXIT
C           IFAIL=1 IF NDIM.LT.1 OR NDIM.GT.20 OR LWA TOO SMALL
C     **************************************************************
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='D01FBF')
C     .. Scalar Arguments ..
      INTEGER                          IFAIL, LWA, NDIM
C     .. Array Arguments ..
      DOUBLE PRECISION                 ABSCIS(LWA), WEIGHT(LWA)
      INTEGER                          NPTVEC(NDIM)
C     .. Function Arguments ..
      DOUBLE PRECISION                 FUN
      EXTERNAL                         FUN
C     .. Local Scalars ..
      DOUBLE PRECISION                 GWT, ONE, ZERO
      INTEGER                          IERROR, ISUMA, ISUMB, J, NDIMM,
     *                                 NVSUM
C     .. Local Arrays ..
      DOUBLE PRECISION                 GSUMS(20), Z(20)
      INTEGER                          ICOUNT(20)
      CHARACTER*1                      P01REC(1)
C     .. External Functions ..
      INTEGER                          P01ABF
      EXTERNAL                         P01ABF
C     .. Data statements ..
      DATA                             ZERO, ONE/0.0D0, 1.0D0/
C     .. Executable Statements ..
C       INTIALISATION AND CHECKING OF PARAMETERS
      IF (NDIM.LT.1 .OR. NDIM.GT.20) GO TO 160
      NDIMM = NDIM - 1
      NVSUM = 0
      DO 20 J = 1, NDIM
         NVSUM = NVSUM + NPTVEC(J)
         GSUMS(J) = ZERO
         ICOUNT(J) = 1
   20 CONTINUE
      IF (NVSUM.GT.LWA) GO TO 160
C       COUNT TO NEXT FUNCTION EVALUATION POINT
   40 IF (NDIMM.EQ.0) GO TO 80
      DO 60 J = 1, NDIMM
         IF (ICOUNT(J).LE.NPTVEC(J)) GO TO 100
         ICOUNT(J) = 1
         GSUMS(J+1) = GSUMS(J+1) + GSUMS(J)
         GSUMS(J) = ZERO
         ICOUNT(J+1) = ICOUNT(J+1) + 1
   60 CONTINUE
   80 IF (ICOUNT(NDIM).GT.NPTVEC(NDIM)) GO TO 140
  100 ISUMA = 0
C       DETERMINTE WEIGHT PRODUCT  AND EVALUATE INTEGRAND
      GWT = ONE
      DO 120 J = 1, NDIM
         ISUMB = ISUMA + ICOUNT(J)
         Z(J) = ABSCIS(ISUMB)
         GWT = GWT*WEIGHT(ISUMB)
         ISUMA = ISUMA + NPTVEC(J)
  120 CONTINUE
      GSUMS(1) = GSUMS(1) + GWT*FUN(NDIM,Z)
      ICOUNT(1) = ICOUNT(1) + 1
      GO TO 40
  140 D01FBF = GSUMS(NDIM)
      IERROR = 0
      GO TO 180
  160 IERROR = 1
  180 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
