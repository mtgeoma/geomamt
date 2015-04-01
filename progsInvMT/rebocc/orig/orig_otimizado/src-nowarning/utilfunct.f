
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     UTILITY FUNCTION
C     - CopyVectorR8() - Copy vector vx to vy  (Real*8)
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE CopyVectorR8(x_1,x_2,vx,y_1,y_2,vy)
      INTEGER x_1,x_2,y_1,y_2
      REAL*8  vx(*),vy(*)

      INTEGER ix,iy
 
      IF (y_2-y_1.NE.x_2-x_1) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR COPY VECTOR !!!'
        STOP
      ENDIF
  
      ix = x_1
      DO iy = y_1,y_2
        vy(iy) = vx(ix)
        ix     = ix + 1
      ENDDO

      RETURN
      END ! CopyVectorR8
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE CopyVectorC16(x_1,x_2,vx,y_1,y_2,vy)
      INTEGER x_1,x_2,y_1,y_2
      COMPLEX*16  vx(*),vy(*)

      INTEGER ix,iy
 
      IF (y_2-y_1.NE.x_2-x_1) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR COPY VECTOR !!!'
        STOP
      ENDIF
  
      ix = x_1
      DO iy = y_1,y_2
        vy(iy) = vx(ix)
        ix     = ix + 1
      ENDDO

      RETURN
      END ! CopyVectorC16
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
   
      SUBROUTINE CopyMatrixR8(x00,x01,x10,x11,nx1,nx2,mx,
     >                        y00,y01,y10,y11,ny1,ny2,my)
      INTEGER x00,x01,x10,x11,y00,y01,y10,y11,nx1,nx2,ny1,ny2
      REAL*8  mx(nx1,nx2),my(ny1,ny2)

      INTEGER ix1,ix2,iy1,iy2

      IF ((y01-y00.NE.x01-x00).or.(y11-y10.NE.x11-x10)) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR COPY MATRIX !!!'
        STOP
      ENDIF

      ix1 = x00
      DO iy1 = y00,y01
        ix2 = x10
        DO iy2 = y10,y11
          my(iy1,iy2) = mx(ix1,ix2)
          ix2 = ix2 + 1
        ENDDO
        ix1 = ix1 + 1
      ENDDO

      RETURN
      END ! CopyMatrixR8
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      SUBROUTINE CopyMatrixI4(x00,x01,x10,x11,nx1,nx2,mx,
     >                        y00,y01,y10,y11,ny1,ny2,my)
      INTEGER x00,x01,x10,x11,y00,y01,y10,y11,nx1,nx2,ny1,ny2
      INTEGER  mx(nx1,nx2),my(ny1,ny2)

      INTEGER ix1,ix2,iy1,iy2

      IF ((y01-y00.NE.x01-x00).or.(y11-y10.NE.x11-x10)) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR COPY MATRIX !!!'
        STOP
      ENDIF

      ix1 = x00
      DO iy1 = y00,y01
        ix2 = x10
        DO iy2 = y10,y11
          my(iy1,iy2) = mx(ix1,ix2)
          ix2 = ix2 + 1
        ENDDO
        ix1 = ix1 + 1
      ENDDO

      RETURN
      END ! CopyMatrixI4
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      SUBROUTINE TransMatrixToVectorR8(mx,np1,np2,s1,s2,e1,e2,vx,sa,ea)
      INTEGER np1,np2,s1,s2,e1,e2,sa,ea
      REAL*8  mx(np1,np2),vx(*)

      INTEGER i_1,i_2,ia

      IF ((e2-s2+1)*(e1-s1+1).NE.(ea-sa+1)) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR COPY MATRIX !!!'
        STOP
      ENDIF
  
      ia = sa
      DO i_1 = s1,e1
        DO i_2 = s2,e2
          vx(ia) = mx(i_1,i_2)
          ia = ia + 1
        ENDDO 
      ENDDO

      END ! 
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
 
      SUBROUTINE ConstantMatrixI4(mx,np1,np2,n1,n2,const_val)
      INTEGER np1,np2,n1,n2
      INTEGER mx(np1,np2),const_val

      INTEGER i_1,i_2
      
      DO i_1 = 1,n1
        DO i_2 = 1,n2
          mx(i_1,i_2) = const_val
        ENDDO
      ENDDO

      RETURN
      END ! ConstantMatrixI4

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
      SUBROUTINE ConstantMatrixR8(mx,np1,np2,n1,n2,const_val)
      INTEGER np1,np2,n1,n2
      REAL*8  mx(np1,np2),const_val

      INTEGER i_1,i_2
      
      DO i_1 = 1,n1
        DO i_2 = 1,n2
          mx(i_1,i_2) = const_val
        ENDDO
      ENDDO

      RETURN
      END ! ConstantMatrixR8

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      SUBROUTINE ConstantMatrixC16(mx,np1,np2,n1,n2,const_val)
      INTEGER np1,np2,n1,n2
      REAL*8  const_val
      COMPLEX*16 mx(np1,np2)

      INTEGER i_1,i_2
      
      DO i_1 = 1,n1
        DO i_2 = 1,n2
          mx(i_1,i_2) = const_val
        ENDDO
      ENDDO

      RETURN
      END ! ConstantMatrixC16

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      SUBROUTINE ConstantVectorI4(ix,n,const_ival)
      INTEGER n,ix(*),const_ival

      INTEGER i
      
      DO i = 1,n
        ix(i) = const_ival
      ENDDO

      RETURN
      END ! ConstantVectorI4

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      SUBROUTINE ConstantVectorR8(vx,n,const_val)
      INTEGER n
      REAL*8  vx(*),const_val

      INTEGER i
      
      DO i = 1,n
        vx(i) = const_val
      ENDDO

      RETURN
      END ! ConstantVectorR8

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
      SUBROUTINE ConstantVectorC16(vx,n,const_val)
      INTEGER n
      REAL*8  const_val
      COMPLEX*16  vx(*)

      INTEGER i
      
      DO i = 1,n
        vx(i) = const_val
      ENDDO

      RETURN
      END ! ConstantVectorC16


C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
 
      SUBROUTINE CumulativeDistance(nx,dx,xdis)
      INCLUDE 'constant.h'
      REAL*8  dx(*),xdis(*)
      INTEGER nx,ix
      
      xdis(1) = D0
      DO  ix = 2,nx+1
       xdis(ix) = xdis(ix-1) + dx(ix-1)
      ENDDO 

      RETURN
      END ! CumulativeDistance

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      SUBROUTINE DistanceBetweenBlocks(nx,dx,cx)
      INCLUDE 'constant.h'
      REAL*8  dx(*),cx(*)
      INTEGER nx,ix
      
      DO  ix = 2,nx
       cx(ix) = dx(ix) + dx(ix-1)
      ENDDO 
      cx(1)    = D2*dx(1)
      cx(nx+1) = D2*dx(nx)

      RETURN
      END ! CumulativeDistance

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      INTEGER FUNCTION SumIndx(Indx,np1,np2,np3,np4,
     >                         s1,s2,s3,s4,e1,e2,e3,e4)
      INTEGER Indx(np1,np2,np3,np4),np1,np2,np3,np4
      INTEGER s1,s2,s3,s4,e1,e2,e3,e4

      INTEGER sumx,i_1,i_2,i_3,i_4
       
      sumx = 0
      SumIndx = 0
      DO i_1 = s1,e1
        DO i_2 = s2,e2
          DO i_3 = s3,e3
            DO i_4 = s4,e4
              sumx = sumx + Indx(i_1,i_2,i_3,i_4)
            ENDDO ! i_4
          ENDDO ! i_3
        ENDDO ! i_2
      ENDDO ! i_1
      SumIndx = sumx
     
      RETURN
      END ! SumIndx

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      INTEGER FUNCTION IndexData(Indx,np1,np2,np3,np4,
     >                           iml,irl,ipl,isl,NRes,NPer,NSta)
      INTEGER iml,irl,ipl,isl,np1,np2,np3,np4,NRes(*),NPer(*),NSta(*)
      INTEGER Indx(np1,np2,np3,np4)

      INTEGER idd,im
      INTEGER SumIndx

      idd = 0
      DO im = 1,iml-1
        idd = idd + SumIndx(Indx,np1,np2,np3,np4,
     >                      im,1,1,1,im,NRes(im),NPer(im),NSta(im))
      ENDDO
      idd = idd + SumIndx(Indx,np1,np2,np3,np4,
     >                    iml,1,1,1,iml,irl-1,NPer(iml),NSta(iml))
      idd = idd + SumIndx(Indx,np1,np2,np3,np4,
     >                    iml,irl,1,1,iml,irl,ipl-1,NSta(iml))
      idd = idd + SumIndx(Indx,np1,np2,np3,np4,
     >                    iml,irl,ipl,1,iml,irl,ipl,isl)
      IndexData = idd

      RETURN
      END  ! IndexData()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      INTEGER FUNCTION SumMatrix_I4(mx,np1,np2,s1,s2,e1,e2)
      INTEGER mx(np1,np2),np1,np2,s1,s2,e1,e2

      INTEGER sumx,i_1,i_2

      IF ((s1.LE.0).OR.(s2.LE.0).OR.(e1.LE.0).OR.(e2.LE.0)) THEN
        SumMatrix_I4 = 0
        GOTO 500
      ENDIF

      sumx = 0
      DO i_1 = s1,e1
        DO i_2 = s2,e2
          sumx = sumx + mx(i_1,i_2)
        ENDDO ! i_2
      ENDDO ! i_1
      SumMatrix_I4 = sumx

500   CONTINUE

      RETURN
      END ! SumMatrix_I4

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      REAL*8 FUNCTION SumVector_R8(n,a)
      INCLUDE 'constant.h'
      INTEGER n
      REAL*8 a(*)

      INTEGER i
      REAL*8  sumx

      sumx = D0
      DO i = 1,n
        sumx = sumx + a(i)
      ENDDO 
      SumVector_R8 = sumx

      RETURN
      END ! SumVector_R8
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Find length of the string
 
      SUBROUTINE Lenb(string,length)
      CHARACTER*(*) string
      INTEGER nstr,istr,length
 
      nstr = len(string)
      DO istr=nstr,1,-1
         IF (string(istr:istr).ne.' ') THEN
            length = istr
            RETURN
         ENDIF
      ENDDO
      length = 0

      RETURN
      END !Subroutine lenb()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Initialize string
 
      SUBROUTINE IniStr(string)
      CHARACTER*(*) string
      INTEGER nstr,istr
 
      nstr = len(string)
      DO istr=1,nstr
        string(istr:istr) = ' '
      ENDDO

      RETURN
      END !Subroutine IniStr

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Search nearest period (same station)

      SUBROUTINE SearchPeriod(DatInx,SenInx,np1,np2,np3,np4,
     >                        iml,irl,ipl,isl,npp,ipb,ipa)
      INTEGER iml,irl,ipl,isl,np1,np2,np3,np4,ipb,ipa,npp
      INTEGER DatInx(np1,np2,np3,np4),SenInx(np1,np2,np3,np4)
     
      INTEGER ip
   
      ipb = 0
      ipa = 0

      DO ip = ipl-1,1,-1
        IF ((SenInx(iml,irl,ip,isl).EQ.1).AND.
     >      (DatInx(iml,irl,ip,isl).EQ.1)) THEN
           ipb = ip
           GOTO 50
        ENDIF
      ENDDO
50    CONTINUE

      DO ip = ipl+1,npp
        IF ((SenInx(iml,irl,ip,isl).EQ.1).AND.
     >      (DatInx(iml,irl,ip,isl).EQ.1)) THEN
           ipa = ip
           GOTO 100
        ENDIF
      ENDDO
100   CONTINUE
 
      RETURN
      END  ! SearchPeriod()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Search nearest station (same period)

      SUBROUTINE SearchStation(DatInx,SenInx,np1,np2,np3,np4,
     >                         iml,irl,ipl,isl,nss,isb,isa)
      INTEGER iml,irl,ipl,isl,np1,np2,np3,np4,isb,isa,nss
      INTEGER DatInx(np1,np2,np3,np4),SenInx(np1,np2,np3,np4)
     
      INTEGER is

      isb = 0
      isa = 0
   
      DO is = isl-1,1,-1
        IF ((SenInx(iml,irl,ipl,is).EQ.1).AND.
     >      (DatInx(iml,irl,ipl,is).EQ.1)) THEN
           isb = is
           GOTO 50
        ENDIF
      ENDDO
50    CONTINUE

      DO is = isl+1,nss
        IF ((SenInx(iml,irl,ipl,is).EQ.1).AND.
     >      (DatInx(iml,irl,ipl,is).EQ.1)) THEN
           isa = is
           GOTO 100
        ENDIF
      ENDDO
100   CONTINUE

      RETURN
      END  ! SearchStation()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Search nearest occupiled sites
C     1 2 3
C     4 X 5
C     6 7 8     

      SUBROUTINE SearchBoth(DatInx,SenInx,SknDepth,Period,StaPos,
     >                      np1,np2,np3,np4,
     >                      area,iml,irl,ipl,isl,npp,nss,ipp,iss)
      INCLUDE 'constant.h'
      INTEGER iml,irl,ipl,isl,np1,np2,np3,np4,ipp,iss,area,npp,nss
      INTEGER DatInx(np1,np2,np3,np4),SenInx(np1,np2,np3,np4)
      REAL*8  SknDepth(np1,np3),Period(np1,np3),StaPos(np1,np4)

      INTEGER ip,is,ip1,ip2,is1,is2,pdir,sdir
      REAL*8  min_dist,dist,dss,dpp

C     Upper Left Side
      IF (area.EQ.1) THEN
        ip1  = ipl - 1
        ip2  = 1
        pdir = -1
        is1  = isl - 1
        is2  = 1
        sdir = -1
      ENDIF     

C     Lower Left Side
      IF (area.EQ.6) THEN
        ip1  = ipl + 1
        ip2  = npp
        pdir = 1
        is1  = isl - 1
        is2  = 1
        sdir = -1
      ENDIF     

C     Upper Right Side
      IF (area.EQ.3) THEN
        ip1  = ipl - 1
        ip2  = 1
        pdir = -1
        is1  = isl + 1
        is2  = nss
        sdir = 1
      ENDIF     

C     Lower Right Side
      IF (area.EQ.8) THEN
        ip1  = ipl + 1
        ip2  = npp
        pdir = 1
        is1  = isl + 1
        is2  = nss
        sdir = 1
      ENDIF     

      min_dist = 3.40282347e+38
      ipp      = 0
      iss      = 0
      DO ip = ip1,ip2,pdir
        DO is = is1,is2,sdir
          IF ((SenInx(iml,irl,ip,is).EQ.1).AND.
     >        (DatInx(iml,irl,ip,is).EQ.1)) THEN
            CALL FindDistance(iml,irl,ipl,isl,ip,is,
     >           SknDepth,Period,StaPos,np1,np3,np4,dist,dss,dpp)
            IF (dist.LT.min_dist) THEN
              min_dist = dist
              ipp      = ip
              iss      = is
            ENDIF
          ENDIF 
        ENDDO ! is
      ENDDO ! ip

      IF (min_dist.EQ.3.40282347e+38) THEN
        min_dist = D0
      ENDIF


      RETURN
      END ! SearchBoth()
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE FindDistance(iml,irl,ipl,isl,ipp,iss,
     >                        SknDepth,Period,StaPos,np1,np3,np4,
     >                        dist,ds,dp)
      INCLUDE 'constant.h'
      INTEGER iml,irl,ipl,isl,np1,np3,np4,ipp,iss
      REAL*8  SknDepth(np1,np3),Period(np1,np3),StaPos(np1,np4)
      REAL*8  ds,dp,dist
 
      REAL*8  dss,dsp

      dist = D0
      dss  = DABS(StaPos(iml,isl) - StaPos(iml,iss))
      dsp  = DABS(SknDepth(iml,ipl) - SknDepth(iml,ipp))
      dist = DSQRT(dss**D2 + dsp**D2)
      
      ds = dss
      dp = DABS(DLOG10(Period(iml,ipl)) - DLOG10(Period(iml,ipp)))

      RETURN
      END ! FindDistance

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Find nearest point
C     algorithms...  go 8 directions..  need more than 1 point
C     1 2 3          
C     4 x 5         
C     6 7 8        
C     region 2 & 7 : same station but before and after period X
C     region 4 & 5 : same period  but before and after station X
C     region 1,3,6 & 8 : different period and station from X


      SUBROUTINE FindNearest(iml,irl,ipl,isl,npp,nss,DatInx,SenInx,
     >                       SknDepth,Period,StaPos,np1,np2,np3,np4,
     >                       io,jo,ds,dp)
      INCLUDE 'constant.h'
      INTEGER iml,irl,ipl,isl,np1,np2,np3,np4,npp,nss
      INTEGER DatInx(np1,np2,np3,np4),SenInx(np1,np2,np3,np4)
      REAL*8  SknDepth(np1,np3),Period(np1,np3),StaPos(np1,np4)
      INTEGER io(8),jo(8)
      REAL*8  ds(8),dp(8)
 
      INTEGER ii,i1,i2,i3,i4,i5,i6,i7,i8,j1,j3,j6,j8
      REAL*8  dist

      DO ii = 1,8
        io(ii) = 0
        jo(ii) = 0
        ds(ii) = D0
        dp(ii) = D0
      ENDDO

      IF (SenInx(iml,irl,ipl,isl).EQ.1) THEN
        RETURN
      ENDIF


C     same station : area 2 & 7
      CALL SearchPeriod(DatInx,SenInx,np1,np2,np3,np4,
     >                  iml,irl,ipl,isl,npp,i2,i7)
      IF (i2.NE.0) THEN
        io(2) = i2
        jo(2) = isl
        CALL FindDistance(iml,irl,ipl,isl,i2,isl,
     >       SknDepth,Period,StaPos,np1,np3,np4,dist,ds(2),dp(2))
      ENDIF
      IF (i7.NE.0) THEN
        io(7) = i7
        jo(7) = isl
        CALL FindDistance(iml,irl,ipl,isl,i7,isl,
     >       SknDepth,Period,StaPos,np1,np3,np4,dist,ds(7),dp(7))
      ENDIF


C     same period : area 4 & 5
      CALL SearchStation(DatInx,SenInx,np1,np2,np3,np4,
     >                   iml,irl,ipl,isl,nss,i4,i5)
      IF (i4.NE.0) THEN
        io(4) = ipl
        jo(4) = i4
        CALL FindDistance(iml,irl,ipl,isl,ipl,i4,
     >       SknDepth,Period,StaPos,np1,np3,np4,dist,ds(4),dp(4))
      ENDIF
      IF (i5.NE.0) THEN
        io(5) = ipl
        jo(5) = i5
        CALL FindDistance(iml,irl,ipl,isl,ipl,i5,
     >       SknDepth,Period,StaPos,np1,np3,np4,dist,ds(5),dp(5))
      ENDIF


C     find 1 & 6
      CALL SearchBoth(DatInx,SenInx,SknDepth,Period,StaPos,
     >     np1,np2,np3,np4,1,iml,irl,ipl,isl,npp,nss,i1,j1)
      IF (i1.NE.0) THEN
        io(1) = i1
        jo(1) = j1
        CALL FindDistance(iml,irl,ipl,isl,i1,j1,
     >       SknDepth,Period,StaPos,np1,np3,np4,dist,ds(1),dp(1))
      ENDIF
      CALL SearchBoth(DatInx,SenInx,SknDepth,Period,StaPos,
     >     np1,np2,np3,np4,6,iml,irl,ipl,isl,npp,nss,i6,j6)
      IF (i6.NE.0) THEN
        io(6) = i6
        jo(6) = j6
        CALL FindDistance(iml,irl,ipl,isl,i6,j6,
     >       SknDepth,Period,StaPos,np1,np3,np4,dist,ds(6),dp(6))
      ENDIF

C     Find 3 & 8
      CALL SearchBoth(DatInx,SenInx,SknDepth,Period,StaPos,
     >     np1,np2,np3,np4,3,iml,irl,ipl,isl,npp,nss,i3,j3)
      IF (i3.NE.0) THEN
        io(3) = i3
        jo(3) = j3
        CALL FindDistance(iml,irl,ipl,isl,i3,j3,
     >       SknDepth,Period,StaPos,np1,np3,np4,dist,ds(3),dp(3))
      ENDIF
      CALL SearchBoth(DatInx,SenInx,SknDepth,Period,StaPos,
     >     np1,np2,np3,np4,8,iml,irl,ipl,isl,npp,nss,i8,j8)
      IF (i8.NE.0) THEN
        io(8) = i8
        jo(8) = j8
        CALL FindDistance(iml,irl,ipl,isl,i8,j8,
     >       SknDepth,Period,StaPos,np1,np3,np4,dist,ds(8),dp(8))
      ENDIF

      RETURN
      END ! FindNearest

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      integer function findstr(str1,str2)
      character*(*) str1, str2
c     returns the position of str2 in str1.  Ignores case.
c     returns 0 if str2 not found in str1

      integer i, j, capdif
      logical same

      capdif= ichar('a')-ichar('A')

      do 20 i= 1, len(str1)-len(str2)+1
         do 10 j=1,len(str2)

            same= str1(i+j-1:i+j-1) .eq. str2(j:j)        .or.
    
     &       'A'.le.str2(j:j) .and. str2(j:j).le.'Z' .and.
     &       ichar(str1(i+j-1:i+j-1)) .eq. ichar(str2(j:j))+capdif .or.
 
     &       'a'.le.str2(j:j) .and. str2(j:j).le.'z' .and.
     &       ichar(str1(i+j-1:i+j-1)) .eq. ichar(str2(j:j)) - capdif
 
            if( .not.same ) go to 20
 10      continue
         findstr=i
         return
 20   continue
 
      findstr=0
      return
      end

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      integer function begwrd(string,iwrd)
      integer iwrd
      character*(*) string
c     Returns the index of the first non-blank character in the iwrd'th
c     non-blank word (word are seperated by spaces, tabs or commas).
c     Returns len if iwrd'th word is not found.
      integer i, nword
      logical wasblk
      intrinsic len

      wasblk=.true.
      nword= 0
      do 100 i=1,len(string)
          if( string(i:i).eq.achar(9) .or.
     &        string(i:i).eq.achar(32) .or.
     &        string(i:i).eq.achar(44)    )then

c             /* current character is blank
              wasblk=.true.
          else
              if(wasblk)  nword= nword + 1
              wasblk= .false.
              if(nword.eq.iwrd)then
                  begwrd= i
                  return
              end if
          end if
  100 continue

      begwrd= len(string)
      return
      end

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      integer function endwrd(string,iwrd)
      integer iwrd
      character*(*) string
c     Returns the index of the last non-blank character in the iwrd'th
c     non-blank word (word are seperated by spaces, tabs or commas).
c     Returns len if iwrd'th word is not found.
      integer i, nword
      logical wasblk
      intrinsic len

c     set default output
      endwrd= len(string)

      wasblk=.true.
      nword= 0
      do 100 i=1,len(string)
          if( string(i:i).eq.achar(9) .or.
     &        string(i:i).eq.achar(32) .or.
     &        string(i:i).eq.achar(44)    )then

c             /* current character is blank
              wasblk=.true.
              if(nword.eq.iwrd) exit
          else
              if(wasblk) nword= nword + 1
              wasblk= .false.
              if(nword.eq.iwrd) endwrd= i
          end if
  100 continue

      return
      end

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE ModelChange(Nzb,Ny,RhoOld,RhoNew,RhoChange,
     >           MNormOld,MNormNew,MNormChange)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nzb,Ny
      REAL*8  RhoOld(NZ0MX,NY0MX),RhoNew(NZ0MX,NY0MX),RhoChange
      REAL*8  MNormOld,MNormNew,MNormChange

      INTEGER iz,iy,mm
      REAL*8  rdiff,rdiff2
  
      rdiff  =  D0
      rdiff2 =  D0
      DO iy = 1,Ny
        DO iz = 1,Nzb
          rdiff = DABS(RhoNew(iz,iy)-RhoOld(iz,iy))/RhoOld(iz,iy)
          rdiff2 = rdiff2 + rdiff**D2
        ENDDO ! iz
      ENDDO ! iy

      mm = Nzb*Ny
      RhoChange   = DSQRT(rdiff2/mm)*D100
      MNormChange = ((MNormOld-MNormNew)/MNormNew)*D100
 
      RETURN
      END ! ModelChange
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
