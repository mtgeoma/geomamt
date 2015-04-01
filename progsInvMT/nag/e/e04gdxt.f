      SUBROUTINE E04GDX(IFLAG,M,N,FJAC,LJ,VT,LVT,FVEC,S,WANTVC,UTF,W,LW)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-790 (DEC 1989).
C
C     **************************************************************
C
C     E04GDX (SVDCMP) IS THE NPL NOSL VERSION OF S. HAMMARLING'S
C     ROUTINE F02WAX (SVDGN1).
C
C     .                                          T
C     THE MAIN DIFFERENCES FROM F02WAX ARE THAT P  IS NOW OUTPUT IN
C     VT AND NOT OVERWRITTEN ON FJAC, AND FVEC IS NOT OVERWRITTEN
C     WHEN WANTVC IS .TRUE..
C
C     E04GDX RETURNS PART OF THE SINGULAR VALUE DECOMPOSITION OF THE
C     M*N (M .GE. N) MATRIX FJAC GIVEN BY
C
C     .                   T
C     FJAC  =  Q * (D) * P ,
C     .            (0)
C
C     WHERE Q AND P ARE ORTHOGONAL MATRICES AND D IS AN N*N DIAGONAL
C     MATRIX WITH NON-NEGATIVE DIAGONAL ELEMENTS, THESE BEING THE
C     SINGULAR VALUES OF FJAC.
C
C     .T
C     P  AND THE DIAGONAL ELEMENTS OF D ARE RETURNED IN VT AND S.
C     .                          T
C     IF WANTVC IS .TRUE.  THEN Q FVEC IS ALSO RETURNED IN UTF.
C
C     THE ARRAY W IS USED FOR WORK SPACE.  IT MUST BE OF LENGTH AT
C     LEAST MAX(5*N-4, N + M*N).  LW MUST BE SET TO THE LENGTH OF W.
C
C     SVEN J. HAMMARLING AND BRIAN T. HINDE.
C     D.N.A.C.S., NATIONAL PHYSICAL LABORATORY, ENGLAND
C
C     **************************************************************
C
C     Modified to call new SVD routines, January 1988.
C
C     .. Scalar Arguments ..
      INTEGER           IFLAG, LJ, LVT, LW, M, N
      LOGICAL           WANTVC
C     .. Array Arguments ..
      DOUBLE PRECISION  FJAC(LJ,N), FVEC(M), S(N), UTF(M), VT(LVT,N),
     *                  W(LW)
C     .. Local Scalars ..
      INTEGER           I, NCOLB, NWHY
C     .. Local Arrays ..
      DOUBLE PRECISION  WW(1)
C     .. External Subroutines ..
      EXTERNAL          F01QCF, F01QDF, F02SWF, F02SYF, F02SXF, F06QFF
C     .. Executable Statements ..
      CALL F06QFF('General',M,N,FJAC,LJ,W(N+1),M)
      NWHY = 0
      CALL F01QCF(M,N,W(N+1),M,W,NWHY)
      IF (WANTVC) THEN
         NCOLB = 1
         DO 20 I = 1, M
            UTF(I) = FVEC(I)
   20    CONTINUE
         CALL F01QDF('Transpose','Separate from A',M,N,W(N+1),M,W,1,UTF,
     *               M,WW,NWHY)
      ELSE
         NCOLB = 0
      END IF
      CALL F06QFF('Upper Triangular',N,N,W(N+1),M,VT,LVT)
      CALL F02SWF(N,VT,LVT,S,W,NCOLB,UTF,M,.FALSE.,WW,1,NWHY)
      CALL F02SXF(N,VT,LVT,0,WW,1,W(N+1),NWHY)
      NWHY = IFLAG
      CALL F02SYF(N,S,W,NCOLB,UTF,M,0,WW,1,N,VT,LVT,W(N+1),NWHY)
      IF (NWHY.NE.0) RETURN
      IFLAG = 0
      RETURN
C
C     END OF E04GDX (SVDCMP)
C
      END
