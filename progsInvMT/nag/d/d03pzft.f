      SUBROUTINE D03PZF(NPDE,M,U,NPTS,X,XP,INTPTS,ITYPE,UOUT,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     -----------------------------------------------------
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Routine to compute solution and possibly its first deriv
C     at the interpolation points stored in the array XP.
C     Linear interpolation is used.
C
C     Parameter list :
C
C     XP(INTPTS) ;  the spatial interpolation points,
C                 XP(I) < XP(I+1) , I = 1, IPTS-1
C
C     UOUT(NPDE,IPTS,ITYPE) ;  array to hold the values found by
C                         linear interpolation
C
C           UOUT(J,K,1)  ;  holds U(XP(K),T) for PDE J
C           UOUT(J,K,2)  ;  holds DU/DX of UOUT(J,K,1)
C
C     X(NPTS)    ;  is an array that contains the original mesh
C
C     U(NEQN)    ;  holds the original solution; the first NPDE*NPTS
C                   components of this contains the PDE solution at the
C                   mesh points X(NPTS) and the last NV components are
C                   the additional coupled ode variables
C
C     ITYPE      ;  has value 1 or 2 depending on how many components
C                   of the array UP are required
C
C     IFAIL    ;  error flag ; set to I if XP(I) lies outside the range
C                   (X(1) , X(NPTS)) .
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     ------------------------------------------------------------------
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D03PZF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, INTPTS, ITYPE, M, NPDE, NPTS
C     .. Array Arguments ..
      DOUBLE PRECISION  U(NPDE,NPTS), UOUT(NPDE,INTPTS,ITYPE), X(NPTS),
     *                  XP(INTPTS)
C     .. Scalars in Common ..
      INTEGER           NIJ
C     .. Local Scalars ..
      DOUBLE PRECISION  TWOU
      INTEGER           I, I1, I2, IFAIL1, IFLAG, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(5)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          D03PZW
C     .. Common blocks ..
      COMMON            /AD03PZ/NIJ
C     .. Executable Statements ..
C
      IFLAG = 0
      TWOU = X02AJF()
C
      IF (NPDE.LE.0) THEN
         NREC = 2
         WRITE (P01REC,FMT=99996) NPDE
         IFLAG = 1
         GO TO 60
      END IF
C
      IF (NPTS.LE.2) THEN
         NREC = 2
         WRITE (P01REC,FMT=99995) NPTS
         IFLAG = 1
         GO TO 60
      END IF
C
C
      IF (INTPTS.LE.0) THEN
         NREC = 2
         WRITE (P01REC,FMT=99994) INTPTS
         IFLAG = 1
         GO TO 60
      END IF
C
      IF ((X(NPTS)-X(1)).LT.(TWOU*(NPTS-1.D0))) THEN
         I1 = 1
         I2 = NPTS
         NREC = 5
         IFLAG = 1
         WRITE (P01REC,FMT=99993) I1, X(I1), I2, X(I2)
         GO TO 60
      END IF
C
      DO 20 I = 2, NPTS
         IF ((X(I)-X(I-1)).LT.TWOU) THEN
            I1 = I - 1
            I2 = I
            NREC = 5
            IFLAG = 1
            WRITE (P01REC,FMT=99993) I1, X(I1), I2, X(I2)
            GO TO 60
         END IF
   20 CONTINUE
C
      IF ((XP(INTPTS)-XP(1)).LT.(TWOU*(INTPTS-1.D0))) THEN
         I1 = 1
         I2 = INTPTS
         NREC = 5
         IFLAG = 2
         WRITE (P01REC,FMT=99992) I1, XP(I1), I2, XP(I2)
         GO TO 60
      END IF
C
      DO 40 I = 2, INTPTS
         IF ((XP(I)-XP(I-1)).LT.TWOU) THEN
            I1 = I - 1
            I2 = I
            NREC = 5
            IFLAG = 2
            WRITE (P01REC,FMT=99992) I1, XP(I1), I2, XP(I2)
            GO TO 60
         END IF
   40 CONTINUE
C
      IF (M.LT.0 .OR. M.GT.2) THEN
         NREC = 2
         WRITE (P01REC,FMT=99998) M
         IFLAG = 1
         GO TO 60
      END IF
C
      IF (ITYPE.LT.1 .OR. ITYPE.GT.2) THEN
         NREC = 2
         WRITE (P01REC,FMT=99997) ITYPE
         IFLAG = 1
         GO TO 60
      END IF
C
      IFAIL1 = 0
      CALL D03PZW(XP,UOUT,INTPTS,X,M,U,NPTS,NPDE,ITYPE,IFAIL1)
C
      IF (IFAIL1.EQ.3) THEN
         IFLAG = 3
         NREC = 2
         WRITE (P01REC,FMT=99999) NIJ, XP(NIJ)
         GO TO 60
      END IF
C
   60 CONTINUE
      IFAIL = P01ABF(IFAIL,IFLAG,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, interpolating point',I16,' with the value',
     *       D13.5,/' ** is outside the X range. ')
99998 FORMAT (' ** On entry, M .ne.  0, 1, or 2, ',/' ** M =',I16)
99997 FORMAT (' ** On entry, ITYPE .ne. 1 or 2, ',/' ** ITYPE =',I16)
99996 FORMAT (' ** On entry, NPDE .le. 0, ',/' ** NPDE =',I16)
99995 FORMAT (' ** On entry, NPTS .le. 2, ',/' ** NPTS =',I16)
99994 FORMAT (' ** On entry, INTPTS .le. 0, ',/' ** INTPTS =',I16)
99993 FORMAT (' ** On entry, mesh points are not in ',/' ** strictly i',
     *       'ncreasing order i.e. mesh point no',I16,/' ** with value',
     *       D13.5,/' ** is greater than mesh point',I16,/' ** with va',
     *       'lue',D13.5)
99992 FORMAT (' ** On entry, interpolation points are not in ',/' ** s',
     *       'trictly increasing order i.e. point no',I16,/' ** with v',
     *       'alue',D13.5,/' ** is greater than point',I16,/' ** with ',
     *       'value',D13.5)
      END
