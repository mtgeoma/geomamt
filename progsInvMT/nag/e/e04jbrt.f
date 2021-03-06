      SUBROUTINE E04JBR(N,NFREE,JFIX,GTP,PA,GPJOLD)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     WHEN THE VARIABLE INDEXED BY JFIX IN THE FORMER PERMUTATION
C     OF FREE VARIABLES HAS BEEN FIXED ON A BOUND, E04JBR (MODGTP)
C     MOVES THE ENTRIES WITH INDICES GREATER THAN JFIX ONE POSITION
C     BACKWARDS IN THE VECTORS CONTAINING THE SEARCH DIRECTION AND
C     THE PREVIOUS GRADIENT FOR THE SUBSPACE OF FREE VARIABLES. IT
C     ALSO DEDUCTS FROM THE INNER PRODUCT OF THE OLD PROJECTED
C     GRADIENT AND SEARCH DIRECTION THE PRODUCT OF THE ELEMENTS
C     RELATING TO THE NEW FIXED VARIABLE.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, MARGARET H.
C     WRIGHT AND ENID M. R. LONG, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  GTP
      INTEGER           JFIX, N, NFREE
C     .. Array Arguments ..
      DOUBLE PRECISION  GPJOLD(N), PA(N)
C     .. Local Scalars ..
      INTEGER           I, IPLUS
C     .. Executable Statements ..
      GTP = GTP - GPJOLD(JFIX)*PA(JFIX)
C     MK6 MEND APPLIED 21-7-77
      IF (JFIX.GT.NFREE) RETURN
      DO 20 I = JFIX, NFREE
         IPLUS = I + 1
         PA(I) = PA(IPLUS)
         GPJOLD(I) = GPJOLD(IPLUS)
   20 CONTINUE
      RETURN
C
C     END OF E04JBR (MODGTP)
C
      END
