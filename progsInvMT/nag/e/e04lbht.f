      SUBROUTINE E04LBH(N,INEW,TIGHT,ISTATE,NFREE)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     E04LBH (OFFBND) IS CALLED, WHEN X(INEW) IS RELEASED FROM ITS
C     BOUND, TO UPDATE NFREE AND ISTATE AND TO RESET TIGHT TO
C     .FALSE.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, MARGARET H.
C     WRIGHT AND ENID M. R. LONG, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           INEW, N, NFREE
      LOGICAL           TIGHT
C     .. Array Arguments ..
      INTEGER           ISTATE(N)
C     .. Executable Statements ..
      NFREE = NFREE + 1
      ISTATE(INEW) = NFREE
      TIGHT = .FALSE.
      RETURN
C
C     END OF E04LBH (OFFBND)
C
      END
