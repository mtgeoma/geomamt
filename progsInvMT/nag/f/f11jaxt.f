      SUBROUTINE F11JAX(M,N,IND,IST,IWORK)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Purpose
C     =======
C     F11JAX takes an array IND of M integers lying in the range [1,N]
C     and converts it to a set of linked lists of equal valued indices.
C
C     Arguments
C     =========
C     M      (input) INTEGER
C            On entry, the number of integer values.
C
C     N      (input) INTEGER
C            On entry, the largest value occurring in IND.
C
C     IND    (input/output) INTEGER array, dimension (M)
C            On entry, an array of integer values in the range [1,N].
C            On exit, IND(i) holds the address in the array IND of
C            another element which had the same value on input. If all
C            such elements  have already been addressed IND(i) is set
C            equal to minus the value itself.
C
C     IST    (output) INTEGER array, dimension (N)
C            On exit, IST(i) holds the start address in the array
C            IND of the linked list for the value i.
C
C     IWORK  (workspace) INTEGER array, dimension (N)
C
C     ==================================================================
C
C     .. Scalar Arguments ..
      INTEGER           M, N
C     .. Array Arguments ..
      INTEGER           IND(M), IST(N), IWORK(N)
C     .. Local Scalars ..
      INTEGER           I, INDI, IWI
C     .. External Subroutines ..
      EXTERNAL          F06DBF
C     .. Executable Statements ..
C
C     Initialize IST and IWORK.
C
      CALL F06DBF(N,0,IST,1)
      CALL F06DBF(N,0,IWORK,1)
C
C     Find IST and convert entries in IND to the set of linked lists.
C
      DO 20 I = 1, M
         INDI = IND(I)
         IF (IST(INDI).EQ.0) IST(INDI) = I
         IWI = IWORK(INDI)
         IF (IWI.NE.0) IND(IWI) = I
         IWORK(INDI) = I
   20 CONTINUE
C
C     Set the end of each list to minus the value concerned.
C
      DO 40 I = 1, N
         IWI = IWORK(I)
         IF (IWI.NE.0) IND(IWI) = -I
   40 CONTINUE
C
      RETURN
      END
