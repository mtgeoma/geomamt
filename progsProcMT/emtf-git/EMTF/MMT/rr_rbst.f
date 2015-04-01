c______________________________________________________________________
c
      subroutine rr_rbst(xr,xc,lomit,var,nt,nd,nsta,ih,ie,n1,i1,
     &   xp,itmax_stk,itmax_ln,orient,decl,theta,ngrp,igrp,
     &   lpwout,ioupw,pw,z,sig_s,sig_e,nu,l_ref_proj,bref,period,
     &   pwrmin)

ccc   robust array remote reference subroutine

ccc    ==> XR(nd,nt) is input data (not modified)
ccc    ==> XC(nd,nt) is input "cleaned" data (not modified)
ccc    ==> LOMIT is a logical array; true for data vectors to omit
ccc    ==> VAR(nt) is input estimate of noise variances in each channel
ccc    ==> IH,IE arrays of indices for start point of H/E channel groupings
ccc    ==> I1 is list of channel #s to use for reference (could be all ... of just some)
ccc    ==> N1 is # of channels to use for reference
ccc    ==>  orient, decl channel orientations, local declination
ccc    ==>  theta is rotation direction for output      
ccc    ==>  lpwout = .true. if full array TF is to be output on this
ccc         call to rr_rbst
ccc    ==> ioupw(1) is 0 to NOT output full array TF ... otherwise is
ccc          unit # for output
ccc    ==> ioupw(2) is normal field defn for array TF ... = -1 for avg
ccc    ==>  bref defines reference fields for Pw output file
ccc    <== Z is remote ref local site TF for all stations

ccc    ==> XP is work space

      include 'iosize.inc'
      parameter (npcmax = 4,n2max=ntmx,npc = 2
     &      ,npcs = (npc*(npc+1))/2 )
      parameter (nwork=npcmax*(ntfmax+2*n2max+npcmax+100+1)
     &                 +ntfmax*ntmx)

      complex xr(nd,nt),xc(nd,nt),xp(nd,nt),pctf(2,ntmx)
ccc      complex xr(nd,nt),xc(nd,nt),xp(nd,nt),pctf(2,0:ntmx)
      real var(nt),eval(ntmx),orient(*),decl(*),theta
      complex s(nsmx),se(nsmx),u(ntmx,npcmax),work(nwork),b(ntmx2)
     &  ,a(2,2),c(2),y(ntmx2),pw(*),sig_s(2,2,*),sig_e(nchemx,nchemx,*),
     &   z(nchemx,2,*),bref(ntmx,2),xxinv(npc,npc)
      real sediag(ntmx),period
      integer ih(*),ie(*),i1(n1),i0(ntmx),ioupw(2),
     &    ngrp,igrp(ntmx+1,ngrpmx),nu(ngrpmx)
      logical lrbst,lomit(nd),l_ref_proj,lpwout

ccc   STAGE I :  ESTIMATE POLARIZATION PARAMTERS FOR EACH TIME SEGMENT
      iw2 = npc*nd+1
      ns = (nt*(nt+1))/2

ccc   divide omitted data by 10**10
      nomit = 0
      do i = 1,nd
         if(lomit(i)) then
            nomit = nomit + 1
            do j = 1,nt
               xr(i,j) = xr(i,j)*1.e-10
               xc(i,j) = xc(i,j)*1.e-10
            enddo
         endif
      enddo
             
ccc   normalize "cleaned" data before finding PCs
      do j = 1,n1
         scale = 1./sqrt(var(i1(j)))
         i0(j) = j
         do i = 1,nd
            iw = iw2+(j-1)*nd+i-1
            work(iw) = xc(i,i1(j))*scale
         enddo
      enddo
      
ccc   make SDM
      call rbstk2(work(iw2),n1,nd,s,xp,itmax_stk)
ccc   compute principal components of normalized SDM
      call rsp(s,n1,npcmax,u,eval,0,xp,dum)

      if(l_ref_proj) then
         call refproj(n1,i1,nt,npcmax,u,bref,nd,eval)
      endif

ccc   form NPC linear combinations of data ... estiamtes of polarization
ccc    vectors ...
ccc    (NPC = 2 always in practice)
      call lc_dat(work(iw2),nd,n1,u,n1,npc,i0,work)

ccc   if requested by setting -  option, omit all FCs for which net SNR
ccc   does not meet minimum threshold
      if(pwrmin .gt. 0 ) call minpwr(work,nd,npc,nt,work(iw2),pwrmin
     &           ,nomit)


ccc   STAGE II :  ESTIMATE FULL ARRAY TF, with full residual matrix
ccc        (Huber Weights, pulling individual channels toward predicted)

ccc   work(1)  is begining of PC data series array
ccc   estimate H, E RR TFs for each MT station
ccc      the integers iw# serve to partition the array WORK
ccc      into the different pieces used inside RBSTREG
ccc      WORK should be dimensioned to allow for NPC*(ND+2*N2+NPC+NB+1)
ccc       complex numbers;  NB is block size inside RBSTREG, currently
ccc       set to 100 (bigger than needed, but who cares? ... this is
ccc       tiny compared to ND)
ccc      copy PC data series array; QR in rbstreg will overwrite
      nd2 = npc*nd
      call ccopy(nd2,work,1,work(iw2),1)
         
      iw2 = npc*nd+1       !
      iw3 = iw2+npc*nd     ! start of array bsv
      iw4 = iw3+nt*2       ! start of array xxinv
      iw5 = iw4+npc*npc    ! start of array work (inside rbstreg)
      call rbstreg(work(iw2),xr,xp,nd,npc,nt,pctf,work(iw3),
     &          work(iw5),se,xxinv,itmax_ln,icvg,nomit)
      if(icvg.lt.0) then
         write(*,*) 'in RR_RBST: did not converge'
         write(*,*) 'icvg= ',icvg
      else
         write(*,*) '# of iterations = ',icvg
      endif

ccc   output array TF with error, signal covariance matrices
      nf = nd-nomit

ccc   adjustment of residual covariance ... needs further
ccc    testing/thought ...
ccc      call adj_res(se,pctf,var,nt,u,n1,i1,work(iw2))

      call wrt_pctf(ioupw(1),nt,ns,pctf,se,xxinv,nf,period)

ccc   this minor correction routine still needs to be checked ...
c     if(l_ref_proj) call refproj_cor(nt,u,pctf,sig_s(1,1,nsta),nu)

ccc   Convert PC TF to Plane wave response with average = normal
      ij = 0
      do j = 1,2
         do i = 1,nt
            ij = ij + 1
            y(ij) = pctf(j,i)
         enddo
      enddo

      call mkbih(b,nt,ih,ie,nsta,ioupw(2))
      c(1) = 1.
      c(2) = 0.
      call anfld(y,nt,nt,b,c,a(1,1),pw)
      c(1) = 0.
      c(2) = 1.
      call anfld(y,nt,nt,b,c,a(1,2),pw(nt+1))

ccc   STAGE III  :  estimate TFs for each desired "group" relative
ccc    to the first two channels in the group ... channels in each
ccc    group are specified in array igrp

ccc   In this last stage we do the final step for robust RR (redescending
ccc    influence curve, downeighting (to zero mostly) all channels in a
ccc    group if enough are deemed "bad" ... 
ccc   After computation of TFs, elements needed for the Z**** file
ccc     are computed, and apparent resistivities, phases and error bars
ccc    are computed

ccc   get error variance from robust TF fit
      call diag(se,nt,sediag)

      if(itmax_ln.gt.0) then
         lrbst = .true.
      else
         ndnt = nd*nt
         call ccopy(ndnt,xr,1,xp,1)
         lrbst = .false.
      endif
ccc   start of work array inside rr_wt
      iw3 = iw2 + nd/2+1
      iref_sta = ioupw(2)
      call rr_wt(lrbst,work,xr,xp,lomit,nt,nd,ngrp,igrp,orient,decl,
     &   theta,work(iw3),sediag,work(iw2),z,sig_s,sig_e,
     &   pctf,nu)

ccc   restore omitted data
      do i = 1,nd
         if(lomit(i)) then
            do j = 1,nt
               xr(i,j) = xr(i,j)*1.e+10
               xc(i,j) = xc(i,j)*1.e+10
            enddo
         endif
      enddo

      return
      end
c______________________________________________________________________
c
      subroutine rr_wt(lrbst,xc,xr,xp,lomit,nt,nd,ngrp,igrp,orient,decl
     &   ,theta,work,var,wts,z,sig_s,sig_e,pctf,nu)

      include 'nstamx.inc'
      integer nt,nd,igrp(ntmx+1,*),nu(*)
      real orient(2,*),decl(*),theta,wts(nd),var(nt)
      complex xc(nd,2),xr(nd,nt),xp(nd,nt),work(*),pctf(2,*),
     &        z(nchemx,2,*),sig_s(2,2,*),
     &        sig_e(nchemx,nchemx,*),sig_r(2,2)
      logical lrbst,lomit(nd)

ccc   ==> input: prediciton channels (PC FCs) in xc
ccc              raw data in xr
ccc              predicted data from Huber weights TF in xp
ccc              estimates of local noise variances in var
ccc              the usual array info stuff (nsta, ih, orient, etc)

ccc   - do final weighted RR estimate
ccc   - change coordinates
ccc   - calculate quantities needed for error bars


ccc   do weights etc. separtely for each station/group
      do k = 1,ngrp
ccc      first calculate weights for each FC vector

         nch = igrp(1,k)
         nche = nch-2
         iw1 = 2*nd+1
         iw2 = iw1 + nch*nd
         iw3 = iw2 + nch*nd
         iw4 = iw3 + 2*nt
         iw5 = iw4 + 4
         iw6 = iw5 + nch*nch

         if(lrbst) then
            call mk_wts(xr,xp,var,lomit,nch,nd,igrp(2,k),nu(k),wts)
ccc         (note still not accounting for "isolated point cleaning"
ccc         in error calculation)
         else
ccc         weights are just ones ...
            do i = 1,nd
               wts(i) = 1.
            enddo
            nu(k) = nd
            nomit = 0
         endif             

ccc      multiply all channels in group (plus predictors by weights)
ccc      (also sorts out desired predicted channels ... so call even
ccc        for non-robust case)
         call mult_wts(xc,work,xp,work(iw1),igrp(2,k),nch,nd,wts)

ccc      do standard LS on this modified data, returning TFs
ccc      (ref to local) in array pctf
ccc      calling with max # of iterations set to zero; just LS fit of
ccc      (already) weighted data
         nomit = nd-nu(k)
         call rbstreg(work,work(iw1),work(iw2),nd,2,nch,pctf
     &      ,work(iw3),work(iw6),work(iw5),sig_r,0,icvg,nomit)

ccc      calculate RR TFs; matrices needed for error scalcs
         call mt_tran_err(pctf,work(iw1),nd,nch,nu(k),sig_r,
     &       nche,z(1,1,k),sig_e(1,1,k),sig_s(1,1,k))
      enddo     !  do k = 1,ngrp

      return
      end
c______________________________________________________________________
c
      subroutine mk_wts(xr,xp,var,lomit,nch,nd,igrp,nu,wts)
      integer nch,nd,igrp(nch)
      include 'nstamx.inc'
      complex xr(nd,*),xp(nd,*)
      real wts(nd),var(*),scale(nchmx)
      logical lomit(nd)
ccc   INPUT: raw data, predicted data, variances, for nch
ccc          channels of data (one station);
ccc          lomit = logical array; true to omit data a priori
ccc   OUTPUT: weights, "cleaned up" data (for cases where bad
ccc           data is in only one channel
ccc           The cleaned data to be used is in XP; the weights are in
ccc           WTS
ccc   HOW?  each channel for which complex errors exceed r1 sds (average
ccc         between real and imaginary parts) is replaced by prediction
ccc         If the number of channels replaced exceeds nbadmx,
ccc         the overall weight for the data vector is set to 0.
ccc   
      parameter (r1 = 2.)
      nbadmx = max(1,nint(float(nch)/3.))

      do j = 1,nch
         scale(j) = sqrt(var(igrp(j)))*r1
      enddo

      nu = 0
      do i = 1,nd
         nbad = 0
         do k = 1,nch
            j = igrp(k)
            if(abs(xr(i,j)-xp(i,j)).gt.scale(k)) then
               nbad = nbad + 1
ccc            use prediction instead of data
            else
ccc            use raw data
               xp(i,j) = xr(i,j)
            endif
         enddo
         if((nbad.gt.nbadmx).or.lomit(i)) then
            wts(i) = 0.0
         else
            nu = nu + 1
            wts(i) = 1.0
         endif
      enddo
      return
      end
c______________________________________________________________________
c
      subroutine mult_wts(xc,xcw,xp,xpw,igrp,nch,nd,wts)
      integer nd,nch,igrp(nch)
      complex xc(nd,2),xcw(nd,2),xp(nd,nch),xpw(nd,nch)
      real wts(nd)

      do i = 1,nd
         do j = 1,2
            xcw(i,j) = wts(i)*xc(i,j)
         enddo
         do k = 1,nch
            j = igrp(k)
            xpw(i,k) = wts(i)*xp(i,j)
         enddo
      enddo
      return
      end
c______________________________________________________________________
c
      subroutine mt_tran_err(pctf,xp,nd,nch,nu,rrinv,nche,z,sig_e,
     &   sig_s)
      include 'nstamx.inc'
      complex pctf(2,nch),z(nchemx,2),sig_e(nchemx,nchemx),htf(2,2),det,
     &    sig_s(2,2),xp(nd,nch),rrinv(2,2),cdotc

ccc   updated version of mt_tran, with improved error estimates
ccc   ( ... mt_tran is WRONG)
ccc   ==> takes Remote TFs as input in pctf, rrinv,
ccc       along with reference, data series
ccc   <== Outputs: (1) Z (nche TFs relating final nche channels to
ccc                    first two),
ccc                (2) SIG_E (residual covariance for RR TFs)
ccc                (3) SIG_S (coherent signal power matrix ...
ccc                     [RH]inv [RR] [HR]inv


      det=pctf(1,1)*pctf(2,2)-pctf(1,2)*pctf(2,1)
      if(abs(det).eq.(0.)) go to 60
      htf(1,1)=pctf(2,2)/det
      htf(2,2)=pctf(1,1)/det
      htf(1,2)=-pctf(1,2)/det
      htf(2,1)=-pctf(2,1)/det

      do i = 1,2
         do j = 1,2
         sig_s(i,j) = 0.0
            do k = 1,2
               do l = 1,2
                  sig_s(i,j) = sig_s(i,j)
     &                 + htf(i,k)*rrinv(k,l)*conjg(htf(j,l))
               enddo
            enddo
         enddo
         do j = 1,nche
            z(j,i) = 0.0
            do k = 1,2
               z(j,i) = z(j,i) + htf(i,k)*pctf(k,j+2)
            enddo
         enddo
      enddo

ccc   calculate residual covariance
ccc   first compute residuals to RR fit
      do j = 1,nche
         do i = 1,nd
            do k = 1,2
               xp(i,j+2) = xp(i,j+2) - z(j,k)*xp(i,k)
            enddo
         enddo
      enddo
ccc   now form residual cross-product matrix
      do i = 1,nche
         do j = 1,i
            sig_e(i,j) = cdotc(nd,xp(1,j+2),1,xp(1,i+2),1)/nu
            if(i.ne.j) sig_e(j,i) = conjg(sig_e(i,j))
         enddo
      enddo
      return

60    continue
      write(*,*) 'Determinant is zero!!!!'
      do i = 1,2
         do j = 1,nche
            z(i,j) = (-999.,-999.)
         enddo
         do j = 1,2
            sig_s(i,j) = (0.,0.)
         enddo
      enddo
      return

      end
c______________________________________________________________________
c
      subroutine wrt_pctf(iou_pctf,nt,ns,pctf,cov,xxinv,nf,period)
      complex pctf(2,nt),cov(ns),xxinv(2,2)
      real period
      integer nf
      write(iou_pctf) period,nf,nt,ns
      write(iou_pctf) pctf,xxinv,cov
      return
      end
ccc_____________________________________________________________________
ccc
      subroutine adj_res(cov,pctf,var,nt,y,n1,i1,a)
ccc   adjusts residual covariance to correct for the effect
ccc   of projection of a portion of noise space into predicted
ccc   data vectors before computation of residuals

ccc   THIS IS AN APPROXIMATE (and ad-hoc) CORRECTION
ccc   Idea:  for the correction assume a diagonal error
ccc    covariance calculate the effect on the residual covariance
ccc    and then add this to the covariance

      complex y(n1,2),cov(*),pctf(2,nt),a(nt,nt)
      real var(nt)
      integer nt,k,l,i1(n1)
      
ccc   construnct array prediction filter
      do j = 1,nt
         do i = 1,nt
            a(i,j) = (0.,0.)
         enddo
      enddo
      do j = 1,n1
         scale = 1./sqrt(var(i1(j)))
         do i = 1,nt
            do l = 1,2
              a(i,j) = a(i,j) + pctf(l,i)*y(i1(j),l)
            enddo
            a(i,j) = a(i,j) * scale
         enddo
      enddo
         
ccc   add correction to error covariance
      ii = 0
      do i = 1,nt
         do j = 1,i
            ii = ii + 1
            cov(ii) = cov(ii) - a(i,j)*var(j) - conjg(a(j,i))*var(j) 
            do l = 1,nt
               cov(ii) = cov(ii) + a(i,l)*var(l)*conjg(a(j,l))
            enddo
         enddo
      enddo 
      return
      end
