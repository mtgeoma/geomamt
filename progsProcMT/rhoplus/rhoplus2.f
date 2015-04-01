                 Program rhoplus

      implicit double precision (a-h, o-z)
      character*17  version
      data version /'Rho+ 2.0 18-05-97'/
c
c  The Rho+ program solves exactly the problem of finding the best-fitting
c     1-D model with MT rho-apparent and impedance phase as input.
c  This version can also bound MT data subject to a specified confidence
c     level providing that the target misfit is above the misfit of
c     the best-fitting model.
c  This version can additionally find offsets in gain (within specified bands
c     relative to other bands) that produce the absolute best misfit.
c
c  The theory with several examples of applications to field data are
c     described in:
c
c     Parker, R.L. and J.R. Booker, Optimal One-Dimensional Inversion and
c       Bounding of Magnetotelluric Apparent Resistivity and Phase Measurements,
c       Phys. Earth Planet. Int, 98, 269-282, 1996.
c
c  Rho+ runs under a command-line interpreter. For a dictionary of command,
c     execute the program and then type "?"
c
c
c      Originally programmed by:
c            Robert L. Parker
c              (with help from Kathy Whaler, Phil Stark and others)
c            IGPP, UCSD, La Jolla CA 92037
c
c      With substantial modifications by:
c            John R. Booker
c            Geophysics, Univ. Of Washington, Seattle WA 98195
c
c      Questions or comments should be addressed to:
c
c            booker@geophys.washington.edu
c      
c ************************************************************************
c     IMPORTANT NOTICE: One feature of this program uses
c     copyrighted subroutines from "Numerical Recipies". If you
c     have a license to use these routines, you may set
c     the flag "nrlic" to .true. Search for ++++++ for a
c     further discussion of this and its impact on the functionality
c     of this program.
c
      logical nrlic
c      data nrlic /.false./
      data nrlic /.true./

c ************************************************************************

c**** NOTICE regarding plot output. This program can generate two
c     types of output for use with plotting routines:
c       matrix: files with generic matrices of results. These can
c             easily be imported into programs like matlab.
c       plot: input scripts for plotting packages. Two are currently
c             implemented plotxy and xyplot. To see how to implement
c****         others, search for CCC@.
c
c$$$$ calls builda, bvls, fndofs, forwrd, getchr, getdat, getone, gtband,
c$$$$       model, plot, polres, range, refine, report, sample, scan
c
      parameter (inmx=200)
      character *80 input(inmx)
      common /store/ input
      common /ndict/ iecho,nin,memory,istate(inmx)
      common /io/ inp,iout,iprint,ierr
c

c If you change mxf or mxxi, you need to search through entire code
c   for other occurences.
      parameter (mxf=101,mxxi=450, mxrw=2*mxf)

      common // a(mxrw*mxxi),d(mxrw),x(mxxi),ma,na
      common /bounds/ bl(mxxi),bu(mxxi)
      common /work/ w(mxxi),zz(mxrw),act(mxrw,mxrw+2),istt(mxxi+1)
      common /shared/ pi,twopi,amu0,radian
      common /observ/ freq(mxf),rho(mxf),drh(mxf),phi(mxf),dph(mxf)
     $  ,erh(mxf),eph(mxf),dwtrh(mxf),dwtph(mxf),sgr(mxf),nfrq
      common /final/ respon(2,mxf),bfresp(2,mxf),bound(2,2*mxf)
      common /discrt/ xi(mxxi),nxi
      common /optiml/ wtmax,fitmin,prob
      common /prolix/ loquor
      common /choice/ a0,state

      character*10 state

      common /bopram/ iband,nbands,nstart(4),nend(4),gain(4)
      logical logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs(4)
      character*3 cecho
      common /flags/ logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs,cecho

      common /save/savrho(mxf),savdrh(mxf),savphi(mxf),savdph(mxf)

      logical frzngn, pltall
      character*75 bugfix
      character*80 root
      character*75 sfile
      character*75 outfil
      character*40 chnull

      data flg1,flg2 /.false.,.true./

c  Read commands
      write(ierr,*)
      write(ierr,*)version
      call scan()

c  check whether commands echo is to be turned off
      call getchr('echo',cecho,need)
      if(need.ge.0)then
        if(need.gt.0)then
          if(cecho.eq.'on')then
            iecho= 1
          elseif(cecho.eq.'off')then
            iecho= -1
          else
            write(ierr,*)'Echo unchanged; bad echo option: '//cecho
          endif
        else
          iecho= -iecho
          if(iecho.lt.0)cecho='off'
          if(iecho.gt.0)cecho='on'
        endif
      endif

c set up amount of output
      call getone('print', verbos, isprnt)
c set default print level
      if(isprnt.le.0)verbos=0
      loquor=verbos

c  Redirect output to root.out if there is a root name
c  (otherwise it goes to stdout (unit 6))

      call getchr('root', root, isroot)
      if (isroot .gt. 0 )then
        length=index(root,' ')-1
        outfil=root(1:length)
        outfil(length+1:length+4)='.out'
        iout=16
        open(unit=iout,file=outfil)
      endif
      write(iout,*)version

      call getchr('debug', chnull, ibug)
      if(ibug.ge.0)then
        loquor=5
        if(cecho.eq.'on')iecho=-1
        call getchr('root', root, isroot)
        if(cecho.eq.'on')iecho=1
        if (isroot .gt. 0 )then
          length=index(root,' ')-1
          bugfix(1:length)=root(1:length)
          bugfix(length+1:length+4)='.dbg'
        else
          bugfix='rplus.dbg'
        endif
        open (unit=9, file=bugfix)
        write(9,*)version
        if(cecho.eq.'on')iecho=1
      endif
c
c  Input the data
c
      call getdat()

c set band offset flags. Get indices of
c  band edges for band gain calculation. See subroutine gtband for
c  explanation of the flags

      call gtband()

c     write(0,*)'offst phjudg usephi phfit',offst,phjudg,usephi,phfit

c  Exclude phase from futher consideration if appropriate
 
         if(.not.phfit)then
            do 35 i=1,nfrq
               dph(i)=0.
               dwtph(i)=0.
 35         continue
         endif

c Apply gain corrections in fixed bands (note gains are log10)
      if(fixgn)then
         do 40 i=nstart(1),nend(1)
            g=10.**gain(1)
            rho(i)=rho(i)*g
            drh(i)=drh(i)*g
 40      continue
         do 50 k=2,nbands
            if(.not.incofs(k))then
               do 45 i=nstart(k),nend(k)
                 g=10.**gain(k)
                 rho(i)=rho(i)*g
                 drh(i)=drh(i)*g
 45            continue
            endif
 50      continue

c Check to see if the gain command has at least one band after
c  the first with a gain of 1.00. If this is not true. All he gains
c  are frozen at their stated levels

         if(nbands.gt.1)then
           frzngn = .not. incofs(2)
           do 52 i=3,nbands
             frzngn=.not.(frzngn .or. incofs(i) )
 52        continue
           if(frzngn)then
             offst=.false.
             write(ierr,*)
             write(ierr,*)'>>> Warning: All gains frozen'
             write(iout,*)
             write(iout,*)'>>> Warning: All gains frozen'
           endif
         endif

         write(iout,*)
         if(logdat)then
           write(iout,*)
     &       'Log10 gains added to Log10(rho) prior to computation:'
           write(iout,'(a,i3,a,i3,a,f7.4,a)')' Band ',nstart(1),' :',
     &                              nend(1),' = ',gain(1)
           do 55 k=2,nbands
             write(iout,'(a,i3,a,i3,a,f7.4)')' Band ',nstart(k),' :',
     &                              nend(k),' = ',gain(k)
 55        continue
         else
           write(iout,*)
     &       'Gains multiplied times rho prior to computation:'
           write(iout,'(a,i3,a,i3,a,f7.4,a)')' Band ',nstart(1),' :',
     &                              nend(1),' = ',10.**gain(1)
           do 58 k=2,nbands
             write(iout,'(a,i3,a,i3,a,f7.4)')' Band ',nstart(k),' :',
     &                              nend(k),' = ',10.**gain(k)
 58        continue
         endif
      endif

      call getchr('surface', state, is)
      if (state(1:3) .eq. 'ins')then
          a0=10.0
      else
          a0=0.0
          state='conducting'
      endif

c  Phase I - find band offsets



c  Save original data and error scales

      if(offst)then
          do 60 i=1,nfrq
            savrho(i)=rho(i)
            savdrh(i)=drh(i)
            savphi(i)=phi(i)
            savdph(i)=dph(i)
 60       continue
      endif

      if(offst.and.nrlic)then
          call fndofs()
      elseif(offst)then
          write(iout,*)' '
          write(iout,*) 'Calculation of optimum band offsets'
     &                  //' disabled because nrlic=.false.'
          write(iout,*)'See comments at start of main program'
          write(ierr,*) 'Calculation of optimum band offsets'
     &                  //' disabled because nrlic=.false.'
          write(ierr,*)'See comments at start of main program'
          phfit=.true.
          offst=.false.
      endif

c  Phase II - find best-fitting models

c  calculate fit to un-modified data

      if(offst)then
c  restore errors in phi if phase are included in final fit
c    and subsequent tasks
         if(phfit)then
            do 70 i=1,nfrq
               dph(i)=savdph(i)
 70         continue
         endif
         do 75 k=2,nbands
            do 75 i=nstart(k),nend(k)
c  restore original data in bands that have calculated offsets
              if(incofs(k))then
                rho(i)=savrho(i)*(10.** (-gain(k)))
                drh(i)=savdrh(i)*(10.** (-gain(k)))
              endif
 75      continue
         fitmin=fit(1)
         if(fixgn)then
          write(iout,'(/a/a,1pg11.4)')
     &     'Report: After applying gains given in ''gain'' command',
     &     '        but before applying optimum gains: Chisq = ', fitmin
         else
          write(iout,'(/a,1pg11.4)')
     &     'Report: Before applying optimum gains: Chisq = ', fitmin
         endif
      else
         fitmin=fit(1)
         write(iout,'(/a,g11.4)')'Report: Chisq:',fitmin
      endif
      call forwrd(nfrq, freq, respon)
      call report()

      if(offst)then
         if(fixgn)then
          write(iout,'(a/a,1pg11.4)')
     &     'Report: After applying gains given in ''gain'' command',
     &     '        but before applying optimum gains: Chisq = ', fitmin
         else
          write(iout,'(/a,1pg11.4)')
     &    'Report: Before applying optimum gains: Chisq = ', fitmin
         endif

c  calculate fit to modified data and final responses

         do 80 k=2,nbands
            do 80 i=nstart(k),nend(k)
               rho(i)=savrho(i)
               drh(i)=savdrh(i)
 80      continue
     
         fitmin=fit(1)
         write(iout,'(/a,1pg11.4)')
     &    '        After applying optimum gains: Chisq = ', fitmin

         call forwrd(nfrq, freq, respon)
         call report()
      endif

      do 100 i=1,nfrq
          bfresp(1,i)=respon(1,i)
          bfresp(2,i)=respon(2,i)
 100  continue
c
c  If matrix form output files requested, write out responses
c
      call getchr('matr', chnull, need)
      if(need.ge.0)then
        if (isroot .gt. 0)then
          length=index(root,' ')-1
          sfile=root(1:length)
          sfile(length+1:length+4)='.rsp'
        else
          sfile='rplus.rsp'
        endif
        open (unit=11, file=sfile)

c       write(11,'(1p,4g12.4,0p,3f9.2)')
        write(11,'(1p,g11.4,g12.4,5g11.4)')
     &      (freq(j),sgr(j)*rho(j),drh(j),bfresp(1,j),
     &      phi(j),dph(j),bfresp(2,j), j=1,nfrq)
        write(iout,'(/a,a)')
     &    'Response summary matrix written to ',sfile(1:lnblnk(sfile))
        close(11)
      endif
c
c  See if D+ model should be calculated
c
      call getchr('model', chnull, need)
      if (need .ge. 0) then
        call refine()
        call builda(ma, nfrq, a, d)
        call bvls(0, ma,nxi,a,d, bl,bu, x, w,act,zz,istt, loopA)
        write(iout,'(/a,1pg11.4)')'Final Refinement: Chisq = ', w(1)**2
        call polres()
        call model()
      endif
c
c  Phase III - find lower/upper bounds on phi or rho

      call range(pltall)

c Phase IV - make plot comand file
c    Note that pltall will be .false. if range did not succeed in
c       finding any bounds to plot. In that case only data and reponses
c       will be plotted.

      call plot(pltall)
c
c     write(iout,*)' '
      write(iout,*)'All requested rho+ computations finished'
      write(ierr,*)' '
      stop 'All requested rho+ computations finished'
      end
c_______________________________________________________________________

      blockdata consts
      implicit double precision (a-h, o-z)
      common /shared/ pi,twopi,amu0,radian
      data pi/3.1415926535898d0/, twopi/6.28318508d0/
      data amu0/12.5663706d-7/, radian/0.01745329252d0/
      end
c______________________________________________________________

      blockdata iounit
      parameter (inmx=200)
      common /io/ inp,iout,iprint,ierr
      common /ndict/ iecho,nin,memory,istate(inmx)
      data inp/5/, iout/6/, iprint /0/, iecho/-1/, ierr/0/
      end
c______________________________________________________________

      subroutine refine()
c$$$$ calls sort
c  Increases resolution of solution vector at transition points;
c  Tosses out points of repeated values at 0 or 1
      implicit double precision (a-h, o-z)
      parameter (mxf=101,mxxi=450, mxrw=2*mxf)
      common // a(mxrw*mxxi),d(mxrw),x(mxxi),ma,na
      common /shared/ pi,twopi,amu0,radian
      common /discrt/ xi(mxxi),nxi
      common /prolix/ loquor
      common /optiml/ wtmax,fitmin,prob
      common /io/ inp,iout,iprint,ierr

      logical logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs(4)
      character*3 cecho
      common /flags/ logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs,cecho
c
      if (loquor .ge. 5) write(9,'(a/100(5(1p,g10.3,0pf6.3)/))')
     $'Solution vector',(xi(j),x(j),j=1,nxi)
c
      kxi=0
c  Keep every point with a transition
      do 1100 j=2, nxi
       if (x(j-1) .ne. x(j)) then
         kxi=kxi + 1
         a(kxi)=xi(j-1)
         if (abs(x(j-1)-x(j)) .eq. 1d0) then
           delxi=0.22*(xi(j) - xi(j-1))
           do 1080 kfil=1, 4
             a(kxi+kfil)=a(kxi) + kfil*delxi
 1080      continue
           kxi=kxi + 4
         endif
        endif
 1100 continue
      kxi=kxi+1
      a(kxi)=xi(nxi)
c  Last point is a branch point!  Extend interval
      if (x(nxi) .ne. 0.0) then
        kxi=kxi+1
        a(kxi)=2.0*xi(nxi)
      endif
c
c  Analyse simple transitions
      do 1200 j=2, nxi-1
        if (x(j) .ne. 0d0 .and. x(j) .ne. 1d0) then
          if (x(j-1) .eq. 1d0) then
            kxi=kxi + 1
            a(kxi)=xi(j-1) + x(j)*(xi(j) - xi(j-1))
          elseif (x(j-1) .eq. 0d0) then
            kxi=kxi + 1
            a(kxi)=xi(j-1) + (1.0- x(j))*(xi(j) - xi(j-1))
          endif
        endif
 1200 continue
c
c  Re-order array then copy into xi and optimal xi
      call sort(kxi, a)
      do 1900 j=2, kxi
        xi(j)=a(j)
        if (xi(j-1)/xi(j).gt. 0.999) xi(j-1)=0.5*(xi(j)+xi(j-2))
 1900 continue
      nxi=kxi
c
      if (loquor .ge. 2 .and. flg3)then
        write(iout,'(/a,i5)')'New number of sampling points:', nxi
      endif
      if(flg1)flg3=.false.
      if (loquor .ge. 5) write(9,'(a/(7g11.4))')
     $' Refined   xi ',(xi(j),j=1,nxi)
c
      return
      end
c_______________________________________________________________________

      subroutine polres()
c$$$$ calls sort
c  Finds poles and zeros of response based on the unit step function
c  respresentation held in common // and /discrt/
c  Then constructs a pole-residue representation of MT admittance
c    c(om) = an0 + sum{n=1 to npol} an(n)/(poles(n)+ i*om)
c  Can be used with QD algorithm to find a model in D+
      implicit double precision (a-h, o-z)
      parameter (mxf=101,mxxi=450, mxrw=2*mxf)
      common // a(mxrw*mxxi),d(mxrw),y(mxxi),ma,na
      common /discrt/ yi(mxxi),nyi
      common /prolix/ loquor
      common /choice/ a0, state
      common /xpand/ poles(2*mxf),zeros(2*mxf),an(2*mxf),an0,npol
      common /local/ x(mxxi),xi(mxxi)
      common /io/ inp,iout,iprint,ierr
      logical in

      character*10 state
c
c  True if  0 < x < 1; neat, huh?
      in(xx)=0d0 .lt. xx .and. xx.lt. 1d0
c
c  Copy solution info to local vectors
      do 1001 n=1, nyi
        x(n)=y(n)
        xi(n)=yi(n)
 1001 continue
      nxi=nyi
c  If last point is branch point force a transition to analytic
      if (x(nxi) .gt. 0.0) then
        nxi=nxi + 1
        x(nxi)=0.0
        xi(nxi)=1.05*xi(nxi-1)
        write(ierr,'(/a)') 
     $  '>>> Model is attempting to reverse surface condition'
        write(iout,'(/a)') 
     $  '>>> Model is attempting to reverse surface condition'
      endif
c
      if (loquor .ge. 3) then
        write(iout,'(a)')' ','Transitions on imaginary axis'
        j=1
        write(iout,'(i4,g12.4,f15.9,a)')j,xi(j),x(j), ' '
        do 1100 j=2, nxi
          if (x(j-1).ne.x(j) .or. in(x(j))) then
            write(iout,'(i4,g12.4,f15.9,a)')j,xi(j),x(j), ' *'
          else
            write(iout,'(i4,g12.4,f15.9,a)')j,xi(j),x(j), ' '
          endif
 1100   continue
      endif
c
      npol=0
      nzer=0
      kxi=0
      x1=x(1)
      x(1)=0.0
      x(nxi+1)=0.0
      xi(nxi+1)=1.05*xi(nxi)
      do 1200 j=2, nxi
c  First get the simple poles
        if (x(j-1) .eq. 1d0 .and. x(j) .eq. 0d0) then
           kxi=kxi + 1
           a(kxi)=xi(j-1)
           npol=npol + 1
           poles(npol)=a(kxi)
        elseif (x(j-1).eq.1d0 .and. in(x(j)) .and. x(j+1).eq.0d0)then
           kxi=kxi+1
           a(kxi)=xi(j-1) + x(j)*(xi(j) - xi(j-1))
           npol=npol + 1
           poles(npol)=a(kxi)
        endif
c  Get simple zeros
        if (x(j-1) .eq. 0d0 .and. x(j) .eq. 1d0) then
          kxi=kxi + 1
          a(kxi)=xi(j-1)
          nzer=nzer + 1
          zeros(nzer)=a(kxi)
        elseif (x(j-1).eq.0d0 .and. in(x(j)) .and. x(j+1).eq.1d0)then
          kxi=kxi + 1
          a(kxi)=xi(j-1) + (1.0- x(j))*(xi(j) - xi(j-1))
          nzer=nzer + 1
          zeros(nzer)=a(kxi)
         endif
c  Pole-zero pair
        if (x(j-1).eq.1d0 .and. x(j+1).eq.1d0 .and. in(x(j))) then
           kxi=kxi+2
           a(kxi-1)=xi(j-1) + x(j)*(xi(j) - xi(j-1))/2.0
           a(kxi)=xi(j) + (1.0- x(j))*(xi(j) - xi(j-1))/2.0
           nzer=nzer + 1
           zeros(nzer)=a(kxi)
           npol=npol + 1
           poles(npol)=a(kxi-1)
        endif
c
        if (x(j-1).eq.0 .and. x(j+1).eq.0 .and. in(x(j))) then
           kxi=kxi+2
           a(kxi-1)=xi(j-1) + (1.0- x(j))*(xi(j) - xi(j-1))/2.0
           a(kxi)=xi(j) + x(j)*(xi(j) - xi(j-1))/2.0
           nzer=nzer + 1
           zeros(nzer)=a(kxi-1)
           npol=npol + 1
           poles(npol)=a(kxi)
        endif
c  Nonstandard arrangement
        if (in(x(j-1)) .and. in(x(j))) then
          write(iout,'(/a,g12.4/a)')
     $    '>>> Unusual configuration at xi=',xi(j-1),
     $    '    A pole may be missing in the expansion'
          write(iout,'(4g12.4)') x(j-2),x(j-1),x(j),x(j+1)
        endif
 1200 continue
c
      call sort(kxi, a)
      do 1900 j=1, kxi
        xi(j)=a(j)
 1900 continue
      nxi=kxi
c
c  For surface insulator interchange poles and zeros
      if (a0 .gt. 0d0) then
        do 1950 n=1, npol
          pn=poles(n)
          poles(n)=zeros(n)
          zeros(n)=pn
 1950   continue
        if (x(2) .gt. 0.99999d0) poles(1)=0d0
      else
        if (x(2) .gt. 0d0) zeros(1)=0.0d0
      endif
c
      write(iout,'(/7a)') ('     ',j=1, min(6,npol)) ,'Poles'
      write(iout,'(1p,6g12.5)')(poles(j),j=1, npol)
      write(iout,'(/7a)') ('     ',j=1, min(6,npol)) ,'Zeros'
      write(iout,'(1p,6g12.5)')(zeros(j),j=1, nzer)
c
c  Find Pole-residue expansion
c  Surface conductor first
      bigA=exp(x1/2.0)
      if (a0 .le. 0d0) then
        write(iout,'(//a,1p,g16.5/a)')'   a0 =',a0,
     $'    n    lambda(n)          a(n)'
c  Pole at omega=0
        if (zeros(1) .gt. 0d0) then
          an(1)=bigA
          do 2105 k=1, npol
            an(1)=an(1)*zeros(k)/poles(k)
 2105     continue
            write(iout,'(i5,1p,2g16.5)')0, 0.0,    an(1)
          do 2200 n=1, npol
            an(n+1)=bigA*(poles(n) - zeros(n))/poles(n)
            do 2100 k=1, npol
              if (k .ne. n)
     $      an(n+1)=an(n+1)*(poles(n) - zeros(k))/(poles(n) - poles(k))
 2100       continue
            write(iout,'(i5,1p,2g16.5)')n,poles(n),an(n+1)
 2200     continue
          do 2250 n=npol, 1, -1
            poles(n+1)=poles(n)
 2250     continue
          poles(1)=0d0
          npol=npol + 1
c  No pole omega = 0
        else
          do 2400 n=1, npol
            an(n)=bigA
            do 2300 k=1, npol
              if (k.lt.npol) an(n)=an(n)*(poles(n)-zeros(k+1))
              if (k.ne.n)    an(n)=an(n)/(poles(n)-poles(k))
 2300       continue
            write(iout,'(i5,1p,2g16.5)')n,poles(n),an(n)
 2400     continue
        endif
      else
c  Surface insulator
        a0=bigA
        write(iout,'(//a,1p,g16.5/a)')'   a0 =',a0,
     $'    n    lambda(n)          a(n)'
        do 3400 n=1, npol
          an(n)=-bigA*(poles(n) - zeros(n))
          do 3300 k=1, npol
            if (k.ne.n)  an(n)=an(n)*(poles(n)-zeros(k))/
     $                               (poles(n)-poles(k))
 3300     continue
          write(iout,'(i5,1p,2g16.5)')n,poles(n),an(n)
 3400   continue
      endif
      an0=a0
      return
      end
c_______________________________________________________________________

      subroutine model()
c$$$$ calls qd
      implicit double precision (a-h,o-z)
c  Derives and lists D+ models given partial-fraction
c  representation of admittance in /xpand/
c  Adapted from Dplus system
c
      parameter (mxf=101,mxxi=450, mxrw=2*mxf)
      common /shared/ pi,twopi,amu0,radian
      common /prolix/ loquor
      common /xpand/ poles(2*mxf),zeros(2*mxf),an(2*mxf),an0,npol
c
      common /soln/ al(2,60),nal,nonzer(60)
      common /io/ inp,iout,iprint,ierr
      dimension d(60),tau(60),h(120),c(120)
c
c  Translate variables into local dialect
      emu0=amu0
      nal=npol + 1
      iprint=loquor
      al(1,1)=an0
      al(2,1)=0.0
      do 1010 j=1, npol
        al(1,j+1)=an(j)
        al(2,j+1)=poles(j)
        h(2*j-1)=an(j)
        h(2*j)  =poles(j)
 1010 continue
c
c  Usings Rutishauser's QD algorithm, find equivalent continued fraction
      nal2=2*nal - 2
      call qd(nal-1, h, c(2))
c  Note indexing runs from surface to bottom, contrary to convention
c  of Parker (JGR 1980, vol 85, p4421).  Also applies to model 
c  description.
      c(1)=al(1,1)
      if (al(2,2).eq. 0.0) c(nal2+1)=0.0
      if (iprint .gt. 0) write(iout,'(/a/(1p,8g9.3))')
     $ 'Continued fraction coeffs', c(1),(c(j+1),j=1,nal2)
c  Find and print the actual model from continued fraction coeffs.
      h(1)=c(1)
      d(1)=h(1)
c  Sort out length units:  scale by 1000 for km
      skilo=0.001
      write(iout,'(/a)')'     Units are kilometers and siemens'
      write(iout,'(a)') 'level  depth      conductance   separation'//
     &               '   total conductance'
      totau = 0.

      do 1500 j=2, nal2, 2
        l=j/2
        if (j.gt.2) tau(l)=1.0/(emu0*h(l)*c(j))
        if (j.eq.2) tau(l)=1.0/(emu0*c(j))
        totau = totau +tau(l)
        write(iout, '(i3,4g14.5)') l,skilo*d(l),tau(l),skilo*h(l),totau
        if (c(j+1).eq. 0.0) go to 1500
        h(l+1)=1.0/(emu0*tau(l)*c(j+1))
        d(l+1)=d(l) + h(l+1)
 1500 continue
c  Print out the depth to bottom if this finite
      if (al(2,2).ne.0.0) write(iout,155) skilo*d(l+1),skilo*h(l+1)
 155  format(3x,g14.5,6x,'--  ',3x,2g14.5)
      if (al(2,2).eq.0.0) write(iout,'(a)')
     $' Model terminates with an insulator'
      if (al(2,2).ne.0.0) write(iout,'(a)')
     $' Model terminates with a perfect conductor'
c
c Convert model in D+ to layered conductivity using mid-point
c   distances to adjacent spikes
c   
c  Not implemented
      return
      end
c_______________________________________________________________________

      subroutine qd(n, c, ans)
c$$$$ calls solve
      implicit double precision (a-h,o-z)
      dimension ans(*),c(2,n)
c
c  Program to convert from a partial to a 1/z-type continued fraction
c  based on the Rutishauser QD algorithm
c  in the notation of Bill Gragg's notes:
c         rho(1)    ith rho
c         rho(2)    ith rho primed
c         alamb(1)  (i-1)th lambda
c         alamb(2)  (i-1)th lambda primed
c         alamb(3)  ith lambda
c         alamb(4)  ith lambda primed
c  where i is i of "do 10" loop
c
c  c contains the coefficients of the partial fraction on entry
c  and the coefficients of the continued fraction on exit
c  so the partial fraction coefficients are lost
c
c  alpha array is stored in c(2,i),i=1,n
c  c0 is stored in c(1,1)
c  beta array is stored in c(1,i),i=2,n
c  where n is number of terms in partial fraction
c
c  Internally declared arrays:
c    alamb(4),rho(2)
c
      if(n.eq.1)go to 20
c
      n1=n-1
      do 10 i=2,n
 10   call solve(n, i, c)
c
c  Convert alpha,beta,c0 to rho,lambda,c0 to give fraction in
c  appropriate form. Ans contains ordered c0,rho,lambda arrays
c
 20   ans(1)=c(1,1)
      ans(2)=c(2,1)
      if(n.eq.1)return
      do 30 i=1,n1
        j=2*i+1
        ans(j)=c(1,i+1)/ans(j-1)
        ans(j+1)=c(2,i+1)-ans(j)
 30   continue
      return
      end
c_______________________________________________________________________

      subroutine solve(n, i, c)
c$$$$$ calls update
      implicit double precision (a-h,o-z)
      dimension alamb(4),rho(2),c(2,n)
c
c  Adds in one more partial fraction and computes the new c0
c  and alpha,beta arrays
c
c
c  new values of c0,alpha(1),beta(1)
c
      c0dash=c(1,1)+c(1,i)
      rho(1)=c(2,1)-c(2,i)
      rho(2)=c(1,1)*rho(1)/c0dash
      c(1,1)=c0dash
      alamb(1)=0.d0
      if(i.gt.2)alamb(1)=c(1,2)/rho(1)
      alamb(2)=rho(1)-rho(2)+alamb(1)
      c(2,1)=rho(2)+c(2,i)
      c(1,2)=alamb(2)*rho(2)
c
c  Special case when i=2
c
      if(i.gt.2)go to 10
      c(2,2)=alamb(2)+c(2,i)
      return
c
c  Special case when i=3
c
 10   i1=i-1
      i2=i-2
      if(i.eq.3)go to 20
c
c  Loop through updating successive members of alpha,beta arrays
c
      call update(i2, alamb, c, n, i)
c
c  Calculate final beta value and last two alphas
c
 20   rho(1)=c(2,i1)-c(2,i)-alamb(1)
      rho(2)=alamb(1)*rho(1)/alamb(2)
      alamb(4)=rho(1)-rho(2)
      c(2,i1)=alamb(2)+rho(2)+c(2,i)
      c(1,i)=alamb(4)*rho(2)
      c(2,i)=alamb(4)+c(2,i)
      return
      end
c_______________________________________________________________________

      subroutine update(n, alamb, c, m, i)
c$$$$ calls no other routines
      implicit double precision (a-h,o-z)
      dimension alamb(*),rho(2),c(2,m)
c
c  New values of alpha,beta for old when another partial fraction
c  is added in. Does all but the first and last values.
c  Enter subroutine with alamb(1),alamb(2) set from calculation
c  of alpha(1),beta(1)
c
      do 10 k=2,n
        rho(1)=c(2,k)-c(2,i)-alamb(1)
        rho(2)=alamb(1)*rho(1)/alamb(2)
        alamb(3)=c(1,k+1)/rho(1)
        alamb(4)=rho(1)-rho(2)+alamb(3)
        c(2,k)=alamb(2)+rho(2)+c(2,i)
        c(1,k+1)=rho(2)*alamb(4)
        alamb(1)=alamb(3)
        alamb(2)=alamb(4)
 10   continue
      return
      end
c_______________________________________________________________________

      subroutine sort(n, x)
c$$$$ calls no other routines
c  Re-orders  x  in place to be increasing.
c  Combsort: souped up version of bubble sort 
c  (Lacey & Box, Byte 16, p315, 1991); optimizes better than sort.
      implicit double precision (a-h, o-z)
      dimension x(n)
c
      ngap=n
 1000 ngap=max(int(ngap/1.3), 1)
      if (ngap.eq.9 .or. ngap.eq.10) ngap=11
      isw=0
      do 1500 i=1, n-ngap
        j=i + ngap
        if (x(i) .le. x(j)) goto 1500
        temp=x(i)
        x(i)=x(j)
        x(j)=temp
        isw=1
 1500 continue
      if (isw.eq.1 .or. ngap .gt. 1) goto 1000
      return
      end
c______________________________________________________________________

      function gammln(xx)
c$$$$ calls no other routines
      implicit double precision (a-h, o-z)
c  Finds ln gamma(x) where
c  gamma(x) = int from 0 to inf t**(x-1) * exp(-t) *dt
      dimension cof(6)
      data cof/76.18009173d0, -86.50532033d0, 24.01409822d0,
     $    -1.231739516d0, .120858003d-2, -.536382d-5/
      data stp/2.50662827465d0/
      data half, one, fpf/0.5d0, 1.0d0, 5.5d0/
c
      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*log(tmp)-tmp
      ser=one
      do 1100 j=1, 6
        x=x+one
        ser=ser+cof(j)/x
1100  continue
      gammln=tmp+log(stp*ser)
      return
      end
c_______________________________________________________________________

      function chipt(nu, prob)
c$$$$ calls chiq
      implicit double precision (a-h, o-z)
c  For chi-squared distribution with nu degrees of freedom, finds value
c  of chi-squared that is exceeded with probability prob.
c  Reliable only if prob < 0.5; accuracy 1 in 10**4.
c  Uses Newton's method in the log of Q=1-P
      common /hideq/ q
c
      chisq=nu+sqrt(18.0*nu)
      do 1200 newt=1, 4
        f=log(chiq(nu, chisq)/prob)
        df=-1.0/(q*chisq)
        chisq=chisq - f/df
 1200 continue
      chipt=chisq
      return
      end
c______________________________________________________________________

      function chiq(nu, chisq)
c$$$$ calls gammln
      implicit double precision (a-h, o-z)
c  Complementary chisq probability Q, by continued fraction:
c  Q = 1 - P = chiq = Prob(chi**2 > chisq)
c  = int (chisq, inf) t**(nu/2-1)*exp(-t/2) dt/(2**(nu/2)*gamma(nu/2)
c  See Abramowitz and Stegun Chap 26: 26.4.10
c  Accurate only when chisq > nu.
      common /hideq/ q
c
      a=nu/2.0
      b=chisq/2.0
      if (chisq .lt. nu)
     $print *,'chiq may be unreliable if chisq < nu'
      q=1.0/b
      lmx=max(15.0, sqrt(real(nu)))
      do 1100 n = lmx, 1, -1
        q=1.0/(b+(n-a)/(1.0+n*q))
 1100 continue
      chiq=q*exp(-b+a*log(b)-gammln(a))
      return
      end
c_____________________________________________________________________

      subroutine mult(ma, na, nb, a, b, c)
      implicit double precision (a-h, o-z)
c$$$$ calls no other routines
c  Forms  c  the product of matrices  a , b:  C = A B
      dimension a(ma,na),b(na,nb),c(ma,nb)
c
      do 1200 i=1,ma
        do 1200 j=1,nb
        cij=0.0
        do 1100 k=1,na
          cij=a(i,k)*b(k,j) + cij
 1100   continue
        c(i,j)=cij
 1200 continue
      return
      end
c_____________________________________________________________________

      subroutine span(n, x, pos, x1, x2)
      implicit double precision (a-h, o-z)
c$$$$ calls no other routines
c  Finds minimum x1, and maximum x2, of a vector x for those indices
c  where pos .gt. 0; if all of pos .le. 0 returns first element of x.
      dimension x(*),pos(*)
c
      x1=1.0d39
      x2=-x1
      npos=0
      do 1100 j=1, n
        if (pos(j) .gt. 0d0) then
          npos=npos + 1
          x1=min(x(j), x1)
          x2=max(x(j), x2)
        endif
 1100 continue
      if (npos .eq. 0) then
        x1=x(1)
        x2=x1
      endif
      return
      end
c_______________________________________________________________________

      subroutine range(pltall)

      implicit double precision (a-h, o-z)
      logical pltall
c       is an output flag that is true if bounds are successfully found
c             and is false if no bounds found for any reason
c       This flag is used by the plotting routine to determine what
c             to plot.

c$$$$ calls builda, bvls, getone, sample, task, xnext, report, span
c$$$$ calls chipt, chiq, forwrd
c  Finds upper/lower bounds on rhoa or phi at the list of frequencies
c  defined in the task file and saved in common /chore/
c
      parameter (mxf=101,mxxi=450, mxrw=2*mxf)
      common // a(mxrw*mxxi),d(mxrw),x(mxxi),ma,na
      common /work/ w(mxxi),zz(mxrw),act(mxrw,mxrw+2),istt(mxxi+1)
      common /bounds/ bl(mxxi),bu(mxxi)
      common /observ/ freq(mxf),rho(mxf),drh(mxf),phi(mxf),dph(mxf)
     $  ,erh(mxf),eph(mxf),dwtrh(mxf),dwtph(mxf),sgr(mxf),nfrq
      common /optiml/ wtmax,fitmin,prob
      common /shared/ pi,twopi,amu0,radian
      common /discrt/ xi(mxxi),nxi
      common /chore/ freqs(2*mxf),ntask
      common /final/ respon(2,mxf),bfresp(2,mxf),bound(2,2*mxf)
      logical logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs(4)
      character*3 cecho
      common /flags/ logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs,cecho
      common /io/ inp,iout,iprint,ierr
      common /prolix/ loquor
      common /choice/ a0, state
      parameter (maxit=41)
      dimension gdat(maxit+1),chisq(maxit+1),bestr(2,2*mxf)
      character*7 ul(2),quant(2), note(2,2*mxf)*1, mark(2)*1

      character*10 state

      logical failed, istask, skpdat
      character*70 root
      character*75 bfile
      character*40 chnull

      data mark/'#','*'/, fresh/-1.0d12/
      data ul/' lower ',' upper '/, quant/' phase ',' rho-a '/

c Ascertain whether the measured datum at a frequency or period
c   at which a bound is requested will be used or skipped in
c   the calculation of that bound. The datum will, of course,
c   be used for other bounds.
c   The default is to skip such a datum.
c NOT IMPLEMENTED! Lines of code in this routine associated with
c   this option are incomplete and will not fully achieve object.

C     call getchr('skip',chnull,need)
C     if(need.gt.0)then
C       if(chnull.eq.'on')then
C          skpdat=.true.
C       elseif(chnull.eq.'off')then
C          skpdat=.false.
C       else
C          skpdat=.true.
C       endif
C     else
C       skpdat=.false.
C     endif

      skpdat=.false.
c
c  Get list of frequencies and quantities to be bounded
      call task(istask)
c
c  Save optimal model response at frequencies of bounds in bestr
      if(istask)then
        call forwrd(ntask, freqs, bestr)
        if(loquor.ge.4)then
           write(9,'(/a)')'    freqs  best fit rho-a  best fit phase'
           write(9,'(3g12.4)')(freqs(j),bestr(1,j),bestr(2,j),j=1,ntask)
        endif
      endif
c
c  Copy data frequencies down 1 in arrays and data likewise;
c  the top row will contain the constraint.
c  NOTE: this shifting is done for ALL data, errors and weights
c    to prevent subsequent confusion.
      do 1200 nf=nfrq+1, 2, -1
        freq(nf)=freq(nf-1)
        rho(nf)=rho(nf-1)
        drh(nf)=drh(nf-1)
        erh(nf)=erh(nf-1)
        dwtrh(nf)=dwtrh(nf-1)
        phi(nf)=phi(nf-1)
        dph(nf)=dph(nf-1)
        dwtph(nf)=dwtph(nf-1)
        eph(nf)=eph(nf-1)
 1200 continue
      nfrq=nfrq+1

      if(.not. istask)then
        pltall=.false.
        write(iout,'(a)')' ',
     $  '>>> No bounding tasks'
        if(iout.eq.16)then
          write(ierr,'(a)')' ',
     $    '>>> No bounding tasks'
        endif
        return
      else
        pltall=.true.
      endif
c
c  Default target = 95 percent confidence 
      target=0.95
      call getone('criterion', target, null)
c  Convert confidence to chisq, as necessary
      if (target .lt. 1.d0) target=chipt(ma, 1.d0-target)
      prob=1.d0 - chiq(ma, target)
      write(iout,'(/a,f9.3)')'Just acceptable chi-squared level:',target
      write(iout,'(a,f9.3)') 'Confidence level of this misfit:  ',prob
c
      if (target .le. fitmin) then
       pltall=.false.
       write(iout,'(a)')' ','>>> Target chi-squared less than '
     $                   //'smallest possible'
       write(ierr,'(a)')
     $   'Target chisq below minimum possible. Bounds not computed'
       return
      endif
c
      call span(nfrq, drh, drh, drmin, d2)
      call span(nfrq, dph, dph, dpmin, d2)
      if (dpmin .le. 0.0) dpmin=0.01
c
c  Set up main body of array, repeating top row
      drh(1)=0.0
      dph(1)=1.d0
      call sample()
      call builda(ma, nfrq, a, d)
c
      lsum=0
      itsum=0
      key=0
      nb=0
c
c  Run through lower bounds, then upper: lbub=-1 for lower, +1 upper
      do 4500 lupdn=1, 2
        lbub= 2*lupdn - 3
c
c  Run up through frequencies (f < 0 means bound phase, f >0 rhoa)
c  Phase => type=-1, kind=1; Rhoa => type=+1, kind=2
        delta=fresh
        do 4000 kase=1, ntask
          note(lupdn,kase)=' '
          type=sign(1d0, freqs(kase))
          kind=1.501+type/2.0
          freq(1)=abs(freqs(kase))
c  Initialize 1st constraint with best-fitting response and small error
          phi(1)=bestr(2,kase)
          rho(1)=bestr(1,kase)
          dph(1)=-type*0.01*dpmin
          drh(1)= type*0.01*drmin

c  Set error to cause a datum which is at the same frequency as this
c    task to be skipped for this bound calculation

C         if(skpdat)then
C           write(0,*)'inside if ',skpdat,ntask,nfrq
C           fkase=abs(freqs(kase))
C           do 3410 kk=2,nfrq
C             write(0,*)'inside loop ',kk,skpdat
C             write(0,*)kk,kase,freqs(kase),
C    &          fkase,freq(kk),fkase-freq(kk),' ',fkase.eq.freq(kk)
C             if(fkase.eq.freq(kk))then
C               write(0,*)
C    &'got here',freqs(kase),freqs(kase).gt.0d0,freqs(fkase).lt.0d0
C               if(freqs(kase).gt.0.)then
C                 tmpdrh= drh(kk)
C                 drh(kk)=0.
C                 write(0,*)'rho datum skipped at ',fkase
C               elseif(freqs(kase).lt.0.)then
C                 tmpdph= dph(kk)
C                 dph(kk)=0.
C                 write(0,*)'phi datum skipped at ',fkase
C               else
C                 stop 'Error in range. Should not be here A'
C               endif
C             endif
C3410       continue
C         endif
c
          if (loquor.ge.4 .and. delta.eq.fresh) write(9,'(a)') ' ',
     $    'Fresh start ============================='
c
c  Iteration to match target chi-squared: adjust 1st datum
c  until agreement is reached

          do 3500 iter=1, maxit
c 
            failed=.false.
            if (iter .gt. 1) then
              call builda(ma, 1, a, d)
              call bvls(key, ma,nxi,a,d, bl,bu, x, w,act,zz,istt, loopA)
              chisq(iter) =w(1)**2
              lsum=lsum + loopA
              key=1
c  First time, substitute optimal fit at this frequency
            else
              base=bestr(3-kind, kase)
              gdat(iter)=base
c             if(lbub.eq.-1 .and. kind.eq.2)delta=min(delta,0.5*base)
              delta=min(delta,0.5*base)
              guess=base + delta
              chisq(iter)=fitmin
              loopA=0
            endif
c
c  Choose next datum value
            gdat(iter+1) = xnext(iter, gdat, chisq, target, guess)
c
            if (loquor .ge. 4) then
              write(9,*)
              write(9,'(3a,2i5,2g12.4)')
     $          'Iteration, loopA',quant(kind),'chisq',
     $          iter,loopA,rho(1),chisq(iter)
              write(9,'(a,g11.3,2g12.4,g11.3)')
     $         'gdat(i), delta, guess, gdat(i+1) ',
     $          gdat(iter), delta, guess, gdat(iter+1)
            endif
            
c
c  Transform true datum to generic kind to make misfit a
c  monotone increasing function
            rho(1)=base + lbub*(gdat(iter+1)-base)
            if (rho(1) .le. 0.0) rho(1)=0.01/(1.0 - rho(1))
            phi(1)=rho(1)
c  Quit when target misfit is achieved 
            if (abs(chisq(iter)-target) .lt. .001*target) goto 3700
c
 3500     continue
c
c  Iteration fails to converge

          write(iout,'(/a, 1pg12.4,4a/a,g12.4)')
     $    '>>> Iteration failed at freq',abs(freqs(kase)),
     $    ' while seeking',ul(lupdn),'bound on',quant(kind),
     $    '>>> chisq = ',chisq(maxit)
          write(iout,'(/a/a/(1p,2g12.4))')'   Table of iterates',
     $    '   datum      chisq', (gdat(j),chisq(j),j=1,maxit)
          if(loquor.ge.4)then
            write(9,*)
            write(9,'(/a, 1pg12.4,4a/a,g12.4)')
     $    '>>> Iteration failed at freq',abs(freqs(kase)),
     $    ' while seeking',ul(lupdn),'bound on',quant(kind),
     $    '>>> chisq = ',chisq(maxit)
            write(9,'(/a/a/(1p,2g12.4))')'   Table of iterates',
     $    '   datum      chisq', (gdat(j),chisq(j),j=1,maxit)
            write(9,*)
          endif
          call report()
c
          note(lupdn,kase)=mark(lupdn)
          nb=1
          iter=maxit
          failed=.true.
c
c  Normal exit to this point
 3700     bound(lupdn,kase)=rho(1)
          itsum=itsum + iter
c
          if (loquor .ge.2) write(iout,'(a,1pg12.4,3a,g12.4,a,i3)')
     $   'At freq ',abs(freqs(kase)),ul(lupdn),' bound on',
     $    quant(kind),rho(1), '   Iterations ',iter

          if (loquor .ge. 4) write(9,'(a/(2g12.4))')
     $    '   datum      chisq', (gdat(j),chisq(j),j=1, iter)
c
c  Adjust step for next guess; fresh start when kind changes
c    or last bound failed
          if (freqs(kase)*freqs(kase+1) .lt. 0.0 .or. failed) then
             delta=fresh
          else
             delta=gdat(iter)-gdat(1)
          endif

c  Restore error used to skip measured datum at same freqeuncy as task.

C         if(skpdat)then
C           do 3900 kk=2,nfrq
C             if(fkase.eq.freq(kk))then
C               if(freqs(kase).gt.0.)then
C                 drh(i)=tmpdrh
C                 write(0,*)'rho datum restored at ',fkase
C               elseif(freqs(kase).lt.0.)then
C                 dph(i)=tmpdph
C                 write(0,*)'phi datum restored at ',fkase
C               else
C                 stop 'Error in range. Should not be here B'
C               endif
C             endif
C3900       continue
C         endif
c
 4000   continue
 4500 continue
c
      if (loquor .ge. 1) then
        write(iout,'(/(a,f8.1))')
     $  'Mean number of iterations per bound ',itsum/real(2*ntask),
     $  'Mean number of BVLS loops per call  ',lsum/real(itsum),
     $  'Total number of loop calls          ',real(lsum)
      endif

c
c  if matrix form output files requested, open file to write out bounds

      call getchr('matr', chnull, need)
      if(need.ge.0)then
        if(cecho.eq.'on')iecho=-1
        call getchr('root', root, isroot)
        if(cecho.eq.'on')iecho=1
        if (isroot .gt. 0)then
          length=index(root,' ')-1
          bfile=root(1:length)
          bfile(length+1:length+4)='.bnd'
        else
          bfile='rplus.bnd'
        endif
        open (12, file=bfile)
      endif

c  Print summary tables: Phase bounds then Rhoa bounds
 
      if (freqs(1) .lt. 0.0) then
         write(iout,'(12x,a)')' ',
     $   '          Bounds on Phase (degrees)',
     $   ' ',
     $   '  frequency   period       lower       upper'
         nph=0
         do 5000 k=1, ntask
           if (freqs(k) .gt. 0.0) goto 5100
           nph=nph+1
           write(iout,'(12x,1p,4g12.4,2a)')-freqs(k),-1.0/freqs(k),
     $       bound(1,k),bound(2,k),note(1,k),note(2,k)
 5000    continue
 
         if (nb .gt. 0) write(iout,'(12x,a)')
     $   '   Lower (#) or upper (*) bound not accurate'
 
c if there are ONLY phase bounding tasks and
c if matrix form output files requested, write out phase bounds
c first line has number of bounds to follow, probability, chisq
c then lines of form: freq, lower bound, upper bound
 
         if(need.ge.0)then
           write(12,'(i8,2g12.4)')nph,prob,fitmin
             do 5150 k=1, ntask
               if (freqs(k) .gt. 0.0) goto 5175
               write(12,'(3g12.4)')-freqs(k),bound(1,k),bound(2,k)
 5150        continue
           write(12,'(a,2g12.4)')'       0',prob,target
         endif
         go to 5555
 
 5100    k1=k
 
c  if there are BOTH phase and rho data, write out matrix phase bounds
         if(need.ge.0)then
           write(12,'(i8,2g12.4)')nph,prob,fitmin
             do 5155 k=1, ntask
               if (freqs(k) .gt. 0.0) goto 5175
               write(12,'(3g12.4)')-freqs(k),bound(1,k),bound(2,k)
 5155      continue
         endif
 
 5175    continue
 
      else
         k1=1
      endif

      if(logdat)then
        write(iout,'(12x,a)')' ',
     $   '  Bounds on Log10 Apparent Resistivity (ohm-m'
      else
        write(iout,'(12x,a)')' ',
     $   '    Bounds on Apparent Resistivity (ohm-m)'
      endif
        write(iout,'(12x,a)')' ',
     $   '  frequency   period       lower       upper'
         nrh=0
c        do 5500 k=k1, ntask
         do 5500 k=ntask, k1, -1
           nrh=nrh+1
           if(logdat)then
            write(iout,'(12x,1p,4g12.4,2a)') freqs(k), 1.0/freqs(k),
     $       log10(bound(1,k)),log10(bound(2,k)),note(1,k),note(2,k)
           else
            write(iout,'(12x,1p,4g12.4,2a)') freqs(k), 1.0/freqs(k),
     $       bound(1,k),bound(2,k),note(1,k),note(2,k)
           endif
 5500    continue
c
c  Write note beneath if necessary
      if (nb .gt. 0) write(iout,'(12x,a)')
     $' Lower (#) or upper (*) bound not accurate'
c
c if matrix form output requested, write out rhoa bounds
c first line has number of bounds to follow, probability, Chisq
c then lines of form: freq, lower bound, upper bound

      if(need.ge.0)then
         if(nph.eq.0)write(12,'(a,2g12.4)')'       0',prob,fitmin
         write(12,'(i8,2g12.4)')nrh,prob,target
c          do 5550 k=k1, ntask
           do 5550 k=ntask, k1, -1
             if(logdat)then
               write(12,'(3g12.4)')freqs(k),
     &                        log10(bound(1,k)),log10(bound(2,k))
             else
               write(12,'(3g12.4)')freqs(k),bound(1,k),bound(2,k)
             endif
 5550      continue
      endif

 5555 continue
      if(need.ge.0)write(iout,'(/a,a)')
     $  'Bounds written in matrix form to '//bfile(1:lnblnk(bfile))

      return
      end
c_______________________________________________________________________

        subroutine task(istask)

      implicit double precision (a-h, o-z)
      logical istask
c        is true if there are bounds to be calculated

c$$$$ calls getchr, sort
c  Reads list of frequencies (or periods) at which bounds on rhoa or phi
c  are desired.   This is a file containing lines each with the format
c              f [rho] [phi]
c  where f is the frequency (or period according to choice in data file)
c  [rho] is the literal string rho when rho-apparent is to be treated
c  [phi] is the literal string phi for phase.  One or both of
c  these strings must appear or the frequency will be ignored.
c
c  Loads results into /chore/ with f<0 meaning phi bound, f >0 rho bound
      parameter (mxf=101,mxxi=450, mxrw=2*mxf)
      common /observ/ freq(mxf),rho(mxf),drh(mxf),phi(mxf),dph(mxf)
     $  ,erh(mxf),eph(mxf),dwtrh(mxf),dwtph(mxf),sgr(mxf),nfrq
      common /shared/ pi,twopi,amu0,radian
      common /chore/ freqs(2*mxf),ntask
      common /prolix/ loquor
      common /io/ inp,iout,iprint,ierr
c
      character*75 flist, line
c
      ntask=0
      call getchr('task',  flist, itsk)
      if (itsk .le. 0) then
        istask=.false.
c       write(iout,'(a)')' ',
c    $  'No bounding tasks requested - Rho+ run complete'
c       if(iout.eq.16)then
c         write(ierr,'(a)')' ',
c    $    'No bounding tasks requested - Rho+ run complete'
c       endif
c       call plot(.false.)
c       stop
        return
      else
        istask=.true.
      endif
c
c  Use data file to supply list - "implicit"
      if (flist .eq. 'implicit') then
        do 1100 j=1, nfrq
          if (drh(j) .le. 0) then
            ntask=ntask + 1
            freqs(ntask)=freq(j)
          endif
          if (dph(j) .le. 0) then
            ntask=ntask + 1
            freqs(ntask)=-freq(j)
          endif
 1100   continue
c
      else
c  Read a task file or standard input
c  If data specified by period rather than frequency ,continue to do so
        istr=index(flist, '*')
        if (istr.eq.0)
     $  open (unit=8, file=flist, status='OLD', err=5000)
c
        call getchr('period', line, nper)
        do 1200 j=1, mxf
          if (istr.eq.0)read (8,'(a)',end=1300) line
          if (istr.gt.0)read (inp,'(a)',end=1300) line

c  Translate tabs (ascii 9) to spaces
          if (index(line(1:lnblnk(line)),char(9)).gt.0) then
            do 1150 k=1,lnblnk(line)
              if (line(k:k) .eq. char(9)) line(k:k)=' '
 1150       continue
          endif
c  Translate nulls (ascii 0) to spaces
          if (index(line, char(0)).gt.0) then
            do 1155 k=1,75
              if (line(k:k) .eq. char(0)) line(k:k)=' '
 1155       continue
          endif
c Skip blank lines
          if (line .eq. ' ') goto 1200
c
C If line has 0 or negative "frequency", stop reading
          read(line,*)  frq
          if (frq .le. 0.0) goto  1300

          i1=index(line, 'h')
c  Skip lines with neither rho nor phi (terminate read if frq=0)
          if (i1 .eq. 0) then
            read(line,*)  frq
            if (frq .le. 0.0) goto  1300
            goto 1200
          endif
c
          if (nper .ge. 0) frq=1d0/frq
          if (index(line(i1-1:75), 'ph') .ne. 0) then
            ntask=ntask + 1
            freqs(ntask)=-frq
          endif
          if (index(line(i1-1:75), 'rh') .ne. 0) then
            ntask=ntask + 1
            freqs(ntask)=frq
          endif
 1200   continue
        write(iout,'(/a,i4)')'>>> Program truncates list of bounds'
 1300   continue
      endif
c
      if (ntask .eq. 0) then
        istask=.false.
        return
      endif
c
      call sort(ntask, freqs)
c
      if (loquor .ge. 2) write(iout,'(//a/a/(1p,5g12.4))')
     $'Frequencies at which bounds are required ',
     $'(Positive values for rho-a, negative for phi)',
     $(freqs(j),j=1, ntask)
      return
c
 5000 write(iout,*)'>>> Unable to open file: ',flist
      stop
c
      end
c_______________________________________________________________________

      function xnext(iter, gdat, chisq, target, guess)
c$$$$ calls find0
c  Inverse interpolation for x for chisq(x)=target
c  Interpolates in log(chi**2), expecting monotone increase in misfit
c  Keeps a table in arrays gdat, chisq
c  Argument  guess  contains guess at first new value; negative guess
c  means start afresh.
c
      implicit double precision (a-h, o-z)
      parameter (maxit=41)
      dimension gdat(*),chisq(*),zero(maxit)
      common /prolix/ loquor
      save fitmax,kount,ibase
c
      if(loquor.ge.5)write(9,*)'Output from xnext:'
      n=iter
      if (n .eq. 1) fitmax=0.0
      if (n .eq. 1) kount=0
      fitmax=max(chisq(n), fitmax)
c  Use guess for 2nd step if there is one
      if (n.eq. 1 .and. guess.gt. 0.0) then
         if(loquor.ge.5)write(9,'(a)')'Using guess'
         xnext=guess
c  Blind increases until chi-squared exceeds target
      elseif (fitmax .lt. target) then
         fact=1 + 0.005*2.0**min(5, n)
         if(loquor.ge.5)
     &     write(9,'(a,g12.4)')'Using blind increase by factor= ',fact
         xnext=gdat(n)*fact
      else
c  Bracketed solution - use Aitken's method on log chisq
        if (kount .eq. 0) then
          if(loquor.ge.5)write(9,*)'n-1 should be > 0: ',n-1
          zero(1)=log(chisq(n-1)/target)
          zero(2)=log(chisq(n)/target)
          ibase=n-1
          if(loquor.ge.5)then
            m=2
            write(9,'(a)')'Input to find0:'
            write(9,'(a)')'ibase, m, ibase+m+1'
            write(9,*)ibase,m,ibase+m+1
            write(9,'(a)')'gdat(ibase..ibase+m+1)'
            write(9,'(5g12.4)')(gdat(i),i=ibase,ibase+m+1)
            write(9,'(a)')'zero(1..m)'
            write(9,'(5g12.4)')(zero(i),i=1,m)
            write(9,'(a)')'End input to find0:'
          endif
          call find0(2, gdat(ibase), zero)
          kount=3
        else
          zero(kount)=log(chisq(n)/target)
          if(loquor.ge.5)then
            write(9,'(a)')'Input to find0:'
            write(9,'(a)')'ibase, kount, ibase+kount+1'
            write(9,*)ibase,kount,ibase+kount+1
            write(9,'(a)')'gdat(ibase..ibase+kount+1)'
            write(9,'(5g12.4)')(gdat(i),i=ibase,ibase+kount+1)
            write(9,'(a)')'zero(1..kount)'
            write(9,'(5g12.4)')(zero(i),i=1,kount)
            write(9,'(a)')'End input to find0:'
          endif
          call find0(kount, gdat(ibase), zero)
          kount=kount + 1
        endif
        if(loquor.ge.5)
     &     write(9,'(a,3i5)')'ibase,kount,ibase+kount-1',
     &                       ibase,kount,ibase+kount-1 
        if(loquor.ge.5)write(9,'(a,i5)')
     &                    'ibase+kount-1 should be <=',maxit
        xnext=gdat(ibase+kount-1)
      endif
c
      if(loquor.ge.5)then
          write(9,*)'xnext= ',xnext
          write(9,*)'End of output from xnext'
          write(9,*)
      endif
      return
      end
c_______________________________________________________________________

      subroutine find0(m, x, y)
c $$$$$ calls nothing
c   Finds next approximation for a root t, where f(t)=0, in a sequence 
c   x(1), x(2), ... x(m), given function values y(1)=f(x(1)), etc.
c   Inserts next element, x(m+1), of the sequence, on the basis of
c   the preceeding values.  User calls f with this value, sets it
c   in y(m+1) and calls find0 again if necessary.  Thus user calls
c   find0 with m=2, 3, ... until satisfied.
c   New value always lies within current best approximating interval.
c   Initial interval must bracket the root, or routine halts.
c   Method:  Aitken's interpolation, or bisection if that fails.
      implicit double precision (a-h, o-z)
      save belo,abov
      common /prolix/ loquor
      common /io/ inp,iout,iprint,ierr

      dimension x(m+1),y(m)
c
c   Set up initial system with linear interpolation
      if (m .le. 2) then
        if (y(1)*y(2) .gt. 0.0) then
          write(iout, '(a)') ' ',
     $    '>>> Initial interval must bracket root in find0','Sorry'
          write(iout,'(a,i5,2g12.4)')'m<=2 & y(1)*y(2)>0 ',m,y(1),y(2)
          write(ierr, '(a)') ' ',
     $    '>>> Initial interval must bracket root in find0'
          write(ierr,'(a,i5,2g12.4)')'m<=2 & y(1)*y(2)>0 ',m,y(1),y(2)
          stop 'find0 failed because interval did not bracket root'
        endif

        yi=y(2)
        x(3)=(x(1)*yi-x(2)*y(1))/(yi-y(1))
        belo=min(x(1), x(2))
        abov=max(x(1), x(2))
        if(loquor.ge.5)then
            write(9,'(a)')
     $      'find0: m<=2: m, x(1), y(1), x(2), y(2), x(3)= '
            write(9,'(i5,5g12.4)')
     $       m, x(1), y(j), x(2), y(2), x(3)
        endif

      else

        sgn=(x(2)-x(1))*y(1)*y(m)
        if (sgn .ge. 0.0) belo=x(m)
        if (sgn .lt. 0.0) abov=x(m)
        if(loquor.ge.5)
     $    write(9,'(a,2g12.4)')'find0: abov, belo= ',abov,belo
c   Form Aitken estimate of root from sequence
        xi=x(m)
        yi=y(m)
        xi=(x(1)*yi - xi*y(1))/(yi - y(1))
        do 1200 j=2, m-1
          if (yi .eq. y(j)) then
            x(m+1)=(abov + belo)/2.0
            return
          endif
          xisav=xi
          xi=(x(j+1)*yi-xisav*y(j))/(yi-y(j))
          if(loquor.ge.5)then
            write(9,'(a)')
     $      'find0: loop j, x(j), y(j), xisav, yi, xi= '
            write(9,'(i5,5g12.4)')
     $       j, x(j), y(j), xisav, yi, xi
          endif
 1200   continue
c  Emergency bisection
        if (xi.gt.abov .or. xi.lt.belo)then
          if(loquor.ge.5)
     $      write(9,'(a,g12.4)')
     $      'find0: using emergency bisection: xi>abov or <belo ',xi
          xi=(abov + belo)/2.0
        endif
        if(loquor.ge.5)write(9,'(a,g12.4)')'find0: final xi= ',xi
        x(m+1)=xi

      endif

      return
      end
c_______________________________________________________________________

      integer function lnblnk( string )
      character string*(*)
      integer len
      
      do i = len(string), 1, -1
        if( string(i:i).ne.' ' .and. string(i:i).ne.char(0) ) then
          lnblnk = i
          return
        endif
      enddo
      
      lnblnk = 0
      return
      end

c_______________________________________________________________________

      subroutine builda(m, mfrq, a,d)
c$$$$ calls no other routines
c  Generates the array a for bvls to work on and the weighted data
c  vector d of ln(rhoa/omega*mu0) and phase, weighted inversely with
c  uncertainties.
      implicit double precision (a-h, o-z)
      dimension a(m,*),d(m)
      parameter (mxf=101,mxxi=450)
      common /shared/ pi,twopi,amu0,radian
      common /observ/ freq(mxf),rho(mxf),drh(mxf),phi(mxf),dph(mxf)
     $  ,erh(mxf),eph(mxf),dwtrh(mxf),dwtph(mxf),sgr(mxf),nfrq
      common /discrt/ xi(mxxi),nxi
      common /prolix/ loquor
      common /io/ inp,iout,iprint,ierr
      common /choice/ a0, state

      character*10 state
      character*15 dreal,dimag
c
c  First time through only, print the transformed weighted data
      if (loquor .ge. 4) write(9,'(a)')' ',
     $'          Transformed data values (std err=1)',
     $'        f (Hz)            Real d         Imag d'
c
c  Runs through by frequency; each freq may generate 0, 1, or 2 
c  rows of  a  depending on what observations are to be fit
      i=0
      do 1500 nf=1, mfrq
        ome=twopi*freq(nf)
        om2=ome**2
c  Check if there is a rhoa row for this frequency
        dreal='           --'
        if (drh(nf) .gt. 0.0) then
          wt=1.0/log(1.0+drh(nf)/rho(nf))
          i=i+1
          a(i,1)=wt
c  Function for surface insulator
          if (a0 .gt. 0d0) then
            do 1200 j=2, nxi
              a(i,j)=-wt*log(((om2+xi(j)*xi(j-1))**2 +
     $                om2*(xi(j)-xi(j-1))**2 )/(om2+xi(j)**2)**2)
 1200       continue
            d(i)=wt*log(rho(nf)/(ome*amu0))
c  Function for surface conductor
          else
            do 1210 j=2, nxi
              a(i,j)=wt*log(((om2+xi(j)*xi(j-1))**2 +
     $                om2*(xi(j)-xi(j-1))**2 )/(om2+xi(j)**2)**2)
 1210       continue
            d(i)=wt*log(ome*rho(nf)/amu0)
          endif
          write(dreal,'(f15.2)') d(i)
        endif
c  Check if there is a phase row and if so, create it
        dimag='           --'
        if (dph(nf) .gt. 0.0) then
          wt=1.0/(radian*dph(nf))
          i=i+1
          a(i,1)=0.0
c  Function for surface insulator
          if (a0 .gt. 0d0) then
            do 1300 j=2, nxi
              a(i,j)=-wt*atan(ome*(xi(j)-xi(j-1))/(om2+xi(j)*xi(j-1)))
 1300       continue
            d(i)=wt*(phi(nf)-90.0)*radian
          else
c  Function for surface conductor
            do 1310 j=2, nxi
              a(i,j)=wt*atan(ome*(xi(j)-xi(j-1))/(om2+xi(j)*xi(j-1)))
 1310       continue
            d(i)=wt*phi(nf)*radian
          endif
          write(dimag,'(f15.2)') d(i)
        endif
        if (loquor .ge. 5)
     $        write(9,'(i5,1pg12.4,2a)')nf,freq(nf),dreal,dimag
c
 1500 continue
c
      return
      end
c_______________________________________________________________________

      subroutine sample()
c$$$$ calls getone, span
c  From data frequencies, generates a sampling of the imaginary axis
c  on which to discretize the integral representation 
c
c  Also works out array size for builda: m, n
c  Also sets lower, upper bounds on solutions
c
c  Receives data from /observ/ and loads /discrt/, /bounds/
c
      implicit double precision (a-h, o-z)
      parameter (mxf=101,mxxi=450,mxrw=2*mxf)
      common /shared/ pi,twopi,amu0,radian
      common // a(mxrw*mxxi),d(mxrw),x(mxxi),ma,na
      common /observ/ freq(mxf),rho(mxf),drh(mxf),phi(mxf),dph(mxf)
     $  ,erh(mxf),eph(mxf),dwtrh(mxf),dwtph(mxf),sgr(mxf),nfrq
      common /discrt/ xi(mxxi),nxi
      common /bounds/ bl(mxxi),bu(mxxi)

      logical logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs(4)
      character*3 cecho
      common /flags/ logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs,cecho
      common /io/ inp,iout,iprint,ierr
c
c  Determine range of sampling in lambda
      call span(nfrq, freq, freq, f1, f2)
c
c  Fill sample array xi: new bold formula!
      anxi=max(90, 60+2*nfrq)
*     anxi=200 + 4*nfrq
      call getone('nlambda', anxi, non)
      nxi=min(mxxi-1, int(anxi))
      xl1=log(f1/3.0)
      xl2=log(60.0*f2)
      dlg=(xl2-xl1)/nxi
      xi(1)=0.0
      do 1500 j=2, nxi
        xi(j)=exp(xl1+(j-2)*dlg)
 1500 continue
c
c  Ascertain array sizes for builda
      m=0
      do 1600 j=1, nfrq
        if (drh(j) .gt. 0.0) m=m + 1
        if (dph(j) .gt. 0.0) m=m + 1
 1600 continue
      ma=m
      na=nxi
c
      if(flg3)
     $write(iout,'(/(a,i4))')'Number of samples on imaginary axis:',nxi,
     $                  '                Number of rows in A:',m
c
c  Load solution bounds into their vectors
      bl(1)=-500.0
      bu(1)=+500.0
      do 1800 j=2, mxxi
        bl(j)=0.0
        bu(j)=1.0
 1800 continue
c
      return
      end
c_______________________________________________________________________

      subroutine forwrd(nfq, frq, rspn)
c$$$$ calls no other routines
c  Computes the response of best-fitting model at list of frequencies
c  frq, setting rhoa and phase (degs) in rspn.
c  Receives solution vector  x in common //
      implicit double precision (a-h, o-z)
      parameter (mxf=101,mxxi=450, mxrw=2*mxf)
      dimension frq(nfq),rspn(2,nfq)
      common // a(mxrw*mxxi),d(mxrw),x(mxxi),ma,na
      common /shared/ pi,twopi,amu0,radian
      common /discrt/ xi(mxxi),nxi
      common /choice/ a0, state

      character*10 state
c
c  Runs through by frequency; each freq generates a rhoa and phase
c  response
      do 1500 nf=1, nfq
        ome=twopi*abs(frq(nf))
        om2=ome**2
        resrh=x(1)
        if (a0 .gt. 0d0) then
c  Surface insulator
          do 1200 j=2, nxi
            resrh=resrh -
     $      x(j)*log(((om2+xi(j)*xi(j-1))**2 +
     $                om2*(xi(j)-xi(j-1))**2 )/(om2+xi(j)**2)**2)
 1200     continue
          rspn(1,nf)=(amu0*ome)*exp(resrh)
        else
c  Surface condunctor
          do 1210 j=2, nxi
            resrh=resrh +
     $      x(j)*log(((om2+xi(j)*xi(j-1))**2 +
     $                om2*(xi(j)-xi(j-1))**2 )/(om2+xi(j)**2)**2)
 1210     continue
          rspn(1,nf)=(amu0/ome)*exp(resrh)
        endif
c
        resph=0.0
        do 1300 j=2, nxi
          resph=resph +
     $    x(j)*atan(ome*(xi(j)-xi(j-1))/(om2+xi(j)*xi(j-1)))
 1300   continue
        if (a0 .gt. 0d0) then
c  Surface insulator or conductor
          rspn(2,nf)=90.0 - resph/radian
        else
          rspn(2,nf)=resph/radian
        endif
c
 1500 continue
c
      return
      end
c_______________________________________________________________________

      subroutine getdat()
c$$$$ calls getchr, span
c$$$$ uses function monotn
c  Read magnetotelluric data from a disk file or the terminal.
c  Each line of the data file has the format:
c               f rhoa err1 phi err2 [dwtrh dwtph]
c  MANDATORY DATA:
c      f=frequency in Hz (or period in s if 'period' included)
c
c      rhoa=apparent resistivity in ohm-meters
c      (can be log10(rhoa) if command logdata invoked)
c
c      err1=error in rhoa; datum ignorded if err1<=0.
c      (must be error in log10(rhoa) if data are log10(rhoa))
c
c      phi=phase of Z in degrees; datum ignorded if err2<=0.
c
c      err2=error in phase, or zero if phi is ignored at this f.
c
c      note: if rhoa or errors<=0, the datum is ignored
c
c  OPTIONAL DATA:
c      dwtrh=downweight for rhoa
c      dwtph=downweight for phi
c      if wt=0, datum is excluded; otherwise datum error multiplied
c        by the downweight.
c
c  Reading terminates at EOF or a zero for f.
c  Blank lines in the data are ignored.
c
      implicit double precision (a-h, o-z)
      common /shared/ pi,twopi,amu0,radian
      parameter (mxf=101,mxxi=450, mxrw=2*mxf)
      common /observ/ freq(mxf),rho(mxf),drh(mxf),phi(mxf),dph(mxf)
     $  ,erh(mxf),eph(mxf),dwtrh(mxf),dwtph(mxf),sgr(mxf),nfrq
      character*75 resdat, line*132
      character*40 chnull
      logical logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs(4)
      character*3 cecho
      common /flags/ logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs,cecho
      common /io/ inp,iout,iprint,ierr
      integer monotn
      external monotn
c     character*20 cstart, cskip
c
      call getchr('data', resdat, mand)
      if (mand .le. 0) then
        write(iout,*)'>>> Data file name not provided: Rho+ halts'
        stop
      endif
c
      if (resdat(1:1) .ne. '*')
     $    open (unit=7, file=resdat, status='OLD', err=5000)
c
      call getchr('period', chnull, isprd)

      call getchr('logd', chnull, islog)
      if(islog.ge.0)then
        logdat=.true.
      else
        logdat=.false.
      endif
c

      nrho=0
      do 1100 j=1, mxf-1

c  read lines into character variable, skipping over blank lines
 1010   if (resdat(1:1) .ne. '*') read(7,'(a)',end=1200) line
        if (resdat(1:1) .eq. '*') read(inp,'(a)',end=1200) line
        if (line .eq. ' ') goto 1010

c  check to see if line contains end-of-data flag
        read(line, *, err=5101) freq(j)
        if (freq(j) .eq. 0.0) goto 1200

c  try to read in data with down-weights
        read(line, *, err=997) freq(j),rho(j),drh(j),phi(j),dph(j)
     $     ,dwtrh(j),dwtph(j)
        go to 999

c  if that fails, read in data without downweights
 997    read(line, *, err=5101) freq(j),rho(j),drh(j),phi(j),dph(j)
        dwtrh(j)=1.
        dwtph(j)=1.

 999    continue

c convert log10 data
        if(logdat)then
          rho(j)=10.**(rho(j))
          drh(j)=2.30259*drh(j)*rho(j)
        endif

c  set downweights to zero for any rhoa<=0 data or any errors<=0
        if(rho(j).le.0.0)then
            rho(j)=-rho(j)
            sgr(j)=-1.
            dwtrh(j)=0.0
        else
            sgr(j)=1.
            nrho=nrho + 1
        endif
        if(drh(j).le.0.0)dwtrh(j)=0.0
        if(dph(j).le.0.0)dwtph(j)=0.0

c  save original error scales for printing and plotting
c  (note that the sign of the rho error scale is used to store
c    the sign of the resistivity)
        erh(j)=(abs(drh(j))+1.d-20)*sgr(j)
        eph(j)=abs(dph(j))

c  use downweights to re-scale errors
        drh(j)=abs(drh(j))*dwtrh(j)
        dph(j)=abs(dph(j))*dwtph(j)

c  convert periods to freqs as necessary
        if (isprd .ge. 0) freq(j)=1.0/freq(j)
 1100 continue

c  If you get here there are too many frequencies in file
c  for the dimensioning.
c  Program will execute with limited number of frequencies
      write(iout,'(/a,i5,a)')'>>> Program will only use first ',
     $mxf - 1,' frequencies in data file'
      write(iout,*)'>>> To alter, edit parameter mxf and recompile'
      j=mxf
 1200 nfrq=j - 1
c
c  List data
      write(iout,'(a)')' ','Data read from file '//resdat,
     $'    Freq (Hz) rho(ohm-m)  rhoerr    lg10rho lg10err phi(deg)  err
     $   downweight'

      do 1500 j=1, nfrq
        sgr(j)=erh(j)/abs(erh(j))
        rlj=log10(rho(j))
        erlj=0.43429*erh(j)/rho(j)
        if(sgr(j).gt.0)then
         if(abs(erh(j)).gt. 1.d-20)then
          write(iout,
     $      '(i3,1pg10.3,g12.4,g11.3,0p,f7.3,1x,f6.3,f8.2,f8.2,2f6.2)')
     $      j,freq(j),rho(j),abs(erh(j)),rlj,erlj,phi(j),eph(j)
     $      ,dwtrh(j),dwtph(j)
         else
          write(iout,
     $      '(i3,1pg10.3,g12.4,a11,0p,f7.3,a,f8.2,f8.2,2f6.2)')
     $      j,freq(j),rho(j),'  0.000    ',rlj,'  0.000',phi(j),eph(j)
     $      ,dwtrh(j),dwtph(j)
         endif
        else
         if(abs(erh(j)).gt. 1.d-20)then
          write(iout,
     $      '(i3,1pg10.3,g12.4,g11.3,0p,a,f8.2,f8.2,2f6.2)')
     $      j,freq(j),-rho(j),abs(erh(j)),'  ------------',phi(j),eph(j)
     $      ,dwtrh(j),dwtph(j)
         else
          write(iout,
     $      '(i3,1pg10.3,g12.4,a11,0p,a,f8.2,f8.2,2f6.2)')
     $      j,freq(j),-rho(j),'  0.000    ','  ------------',
     $      phi(j),eph(j),dwtrh(j),dwtph(j)
         endif
        endif
 1500 continue

c test that frequencies are monotonic
      inon=monotn(freq,nfrq)
      if(inon.ne.0)then
        if(isprd.ge.0)then
          write(iout,*)' '
          write(iout,*)'Periods in input not monotonic at n= ',inon
          write(ierr,*)'Periods in input not monotonic at n= ',inon
        else
          write(iout,*)' '
          write(iout,*)'Frequencies in input not monotonic at n= ',inon
          write(ierr,*)'Frequencies in input not monotonic at n= ',inon
        endif
        stop '   Rho+ quits'
      endif

      if (nrho .eq. 0) goto 5200
c
      return
c
c  Error exits
 5000 write(ierr,*)'>>> Cannot open file:',resdat
      stop
 5101 write(ierr,*)'>>> Trouble reading from file ',resdat
      write(ierr,*)'    Error in following line:'
      write(ierr,*)'=> '//line
      stop
 5200 stop 'Must be at least one rho-a in input data. Rho+ quits'
c
      end
c_______________________________________________________________________

      integer function monotn(V,n)
      implicit double precision (a-h, o-z)
      integer i, n
      dimension V(n)
 
c     returns 0 if V is monotone
c     otherwises returns index of first non-monotone element
 
      logical increas, decreas
 
      increas=.true.
      decreas=.true.
      do 100 i=2,n
          decreas= decreas.and. v(i).lt.v(i-1)
          increas= increas.and. v(i).gt.v(i-1)
          if( .not.decreas .and. .not.increas )then
              monotn= i
              return
          end if
  100 continue
      monotn= 0
 
      return
      end
c______________________________________________________________________

      subroutine scan()
c$$$$ calls no other routines
c  Input routine for commands 
c
c  Reads from the standard input until: eof, "read", "execute", "continue"
c    or "quit" command is encountered.
c  "read" command switches input to a named file
c
c  Saves lines in the input store  /store/ for later retrieval by
c  getarr, getone or getchr.
c
c  Prints help upon request
c
      parameter (inmx=200)
      character*80 input(inmx),line
      character*80 infil
      character*80 ans
      common /store/ input
      common /ndict/ iecho,nin,memory,istate(inmx)
      common /io/ inp,iout,iprint,ierr
c
      nin=0
      memory=inmx
      write(ierr,'(a)') ' ',
     $'                    =================',' ',
     $'Enter commands for analysis of MT responses (? for help)'
c
      do 1500 l=nin+1, inmx
        read(inp,'(80a)', end=3000) line
 1301   continue

        if(line(1:4).eq.'read')then
          
c  Switch file that commands are read from
          do 1332 i=5,80
            if(line(i:i).ne.' ' .and. line(i:i).ne. '	'
     &         .and. line(i:i).ne.char(0))then
              istart=i
              go to 1333
            elseif(i.eq.80)then
              stop 'Read command requires a file name; Rho+ quits'
            endif
 1332     continue

 1333     do 1334 i=istart+1,80
            if(line(i:i).eq.' ' .or. line(i:i).eq. '	'
     &         .or. line(i:i).eq.char(0) .or. line(i:i).eq.'%')then
              iend=i-1
              go to 1336
            endif
 1334     continue

 1336     infil=line(istart:iend)
c  Close current file if it is not stdin
          if(inp.ne.5)close(inp)
          inp=15
          open(unit=inp,file=infil,status='OLD',err=988)
          write(ierr,*)
     &   'Further commands read from file: '//infil(1:lnblnk(infil))
          go to 1500
        endif


        if (line(1:4).eq.'cont' .or. line(1:4).eq.'exec')go to 2000

c  List a glossary of keywords
        if (line(1:1) .eq. '?' .or. line(1:4).eq.'help')then
           write(iout,*)' '
           write(iout, '(a/(2x,a))')
     $'Enter at least first four letters of commands',
     $'   from the following list:',' ',
     $'     (M) means mandatory information',
     $'     [info] is user supplied information',
     $'     { } encloses optional text or response',
     $'     NOTE: Command case IS important',
     $' '
           write(iout,'(a)')
     $'GENERAL:',
     $'%                  All text after a % is ignored',
     $'?                  Remind me of the command list again',
     $'                    prints 1st help screen and then resumes',
     $'                    reading commands if not from keyboard',
     $'r{esume}           Exits help (ignored if not in help)',
     $'help               same as ?',
     $'execute            Quit reading commands and begin calculations',
     $'continue           same as execute',
     $'q{uit}             Quit reading commands and exit program',
     $'echo {choice}      Echoing of commands to output is {on}{off}',
     $'                    DEFAULT is off; no {choice} implies on',
     $'blank line         Will be skipped',
     $'invalid keyword    If 1st 4 characters not a valid command,',
     $'                    the line will be skipped',
     $' '
           write(ierr,*)'Enter: y if you want INPUT OPTIONS help'
           write(ierr,*)'       r to resume reading commands'
           write(ierr,*)'       RETURN to skip to next help topic'
           read(inp,'(a)',err=2001)ans
           if(ans(1:1).ne.'y' .and. ans(1:1).ne.'Y')then
             if(ans(1:1).eq.'r')then
                go to 1379
             elseif(ans(1:1).eq.char(32))then
                go to 1371
             else
                write(0,*)'Invalid input; help quits'
                line=ans
                write(ierr,'(a)') ' ',
     $'                    =================',' ',
     $'Enter commands for analysis of MT responses (? for help)'
                go to 1301
             endif
           endif
           write(iout,'(a)')
     $' ',
     $'INPUT:',
     $'read [filename]    Switch reading commands to named file NOW',
     $'                    (this is a one-way command)',
     $'data [filename](M) Input file of rho-a and phase values',
     $'                    Frequencies can inc or decr monotonically',
     $'period             Data & task files use period, not frequency',
     $'logdata            Data, errors and gains are log10(rho)',
     $'                    Output parallels this choice',
     $' '
           write(iout,'(a)')
     $' ',
     $'INPUT FILE FORMATS:',
     $'Datafile           freq rho rho_err phi phi_err {dwtrho dwtphi}',
     $'                     dwt=0., err or rho <= 0, ignore datum',
     $'Taskfile           freq {rho (text)} {phi (text)}',
     $'                     freq=0.0 or eof terminates either file',
     $' '

 1371      write(ierr,*)'Enter: y if you want OUTPUT OPTIONS help'
           write(ierr,*)'       r to resume reading commands'
           write(ierr,*)'       RETURN to skip to next help topic'
           read(inp,'(a)',err=2001)ans
           if(ans(1:1).ne.'y' .and. ans(1:1).ne.'Y')then
             if(ans(1:1).eq.'r')then
                go to 1379
             elseif(ans(1:1).eq.char(32))then
                go to 1374
             else
                write(0,*)'Invalid input; help quits'
                line=ans
                write(ierr,'(a)') ' ',
     $'                    =================',' ',
     $'Enter commands for analysis of MT responses (? for help)'
                go to 1301
             endif
           endif
           write(iout,'(a)')
     $' ',
     $'OUTPUT:',
     $'printlevel [#]     Amount of output to provide [0 to 5]',
     $'                    DEFAULT is 0',
     $'root {name}        Root name for output files',
     $'                    DEFAULT output goes to stdout if no root',
     $'                    DEFAULT root for optional files is rplus'
           write(iout,'(a)')
     $'plot {type}        Generate plot file [root].plt of type:',
     $'                    plotxy: Parker''s plotxy',
     $'                    xyplot: UW/Chave/Schultz version of plotxy',
     $'                    xyplot fan{cy}: fancier version of xyplot',
     $'                    DEFAULT is xyplot'
CCC@ duplicate 2nd previous line and edit text when types are added
CCC@ edit DEFAULT text if changed
           write(iout,'(a)')
     $'matrix             Matrix form output [root].rsp & [root].bnd',
     $'debug              Print debugging info to [root].dbg',
     $' '

 1374      write(ierr,*)'Enter: y if you want COMPUTATION OPTIONS help'
           write(ierr,*)'       r to resume reading commands'
           write(ierr,*)'       RETURN to skip to next help topic'
           read(inp,'(a)',err=2001)ans
           if(ans(1:1).ne.'y' .and. ans(1:1).ne.'Y')then
             if(ans(1:1).eq.'r')then
                go to 1379
             elseif(ans(1:1).eq.char(32))then
                go to 1375
             else
                write(0,*)'Invalid input; help quits'
                write(ierr,'(a)') ' ',
     $'                    =================',' ',
     $'Enter commands for analysis of MT responses (? for help)'
                line=ans
                go to 1301
             endif
           endif
           write(iout,'(a)')
     $' ',
     $'COMPUTATION:',
     $'surface [choice]   con{ducting} (DEFAULT) or ins{ulating}',
     $'model              Compute bestfit delta function model',
     $'task [filename]    Task file specifies where bounds are wanted',
     $'                    Tasks independent and can be in any order',
     $'task implicit      Search for any missing or excluded data in',
     $'                    input and provide bounds for them',
     $'criterion [#]      Tolerence for bound calculations. [#] can be',
     $'                    0 to 1. Default is 0.95 (95%). If [#] > 1,',
     $'                    it is interpreted as the target ChiSq'
           write(iout,'(a)')
     $'bands [2*N#] {-#} {p{runephase}} Indices of 1st and last data',
     $'                    in N (1 to 4) NON-overlapping bands to',
     $'                    set or compute gains for rho',
     $'                    {-#} Gain estimation options:',
     $'                    >= 0, missing or illegal',
     $'                       All gains held fixed (DEFAULT)',
     $'                    -1 Use all data for gain and bounds',
     $'                    -2 Use only rho data for gain and bounds',
     $'                    -3 Use rho data for gain, all for bounds',
     $'                    Gain in FIRST BAND read always held FIXED',
     $'                    Phase data outside bands are excluded if',
     $'                       option is 0 or {p{runephase} is present'
           write(iout,'(a)')
     $'gains [N#]         Gains to be applied to bands defined',
     $'                     by bands command PRIOR to any computation',
     $'                   Non-unity gain (non-zero logdata offset)',
     $'                     freezes gain in band',
     $'                   Unity gain (zero logdata offset), EXCEPT in',
     $'                     first band read, enables gain estimation',
     $'                     if enabled by bands command option',
     $'                   DEFAULT is unity gain (zero logdata offset)',
     $'                     (i.e. estimation enabled) in all bands',
     $'nlamba [#]         Force [#] of samples on imaginary axis',
     $' '
 
 1375 continue
C1375      write(ierr,*)'Enter: y if you want ADVANCED OPTIONS help'
C          write(ierr,*)'       r or RETURN to resume reading commands'
C          read(inp,'(a)',err=2001)ans
C          if(ans(1:1).ne.'y' .and. ans(1:1).ne.'Y')then
C            if(ans(1:1).eq.'r')then
C               go to 1379
C            elseif(ans(1:1).eq.char(32))then
C               go to 1379
C            else
C               write(ierr,*)'Invalid input; help quits'
C               line=ans
C               go to 1301
C            endif
C          endif
C          write(iout,'(a)')

C    $' ',
C    $'ADVANCED:',
c Following is NOT IMPLEMENTED (see subroutine range)
c    $'skipdatum {choice} Ignoring misfit to measurement at the same',
c    $'                    frequency as a bounding task is {on}{off}',
c    $'                    DEFAULT is off; no {choice} implies on',

1379    continue
        write(ierr,'(a)') ' ',
     $'                    =================',' ',
     $'Enter commands for analysis of MT responses (? for help)'

        endif

        if (line(1:1) .ne. ' ') then
          if (line(1:4) .eq. 'quit' .or. line(1:1) .eq. 'q') goto 3000

c  Translate tabs to spaces
          j1=index(line, '	')
          if (j1 .gt. 0) then
            do 1400 j=5,80
              if (line(j:j) .eq. '	') line(j:j)=' '
 1400       continue
          endif

c  Translate nulls (ascii 0) to spaces
          j1=index(line, char(0))
          if (j1 .gt. 0) then
            do 1450 j=j1,80
              if (line(j:j) .eq. char(0)) line(j:j)=' '
 1450       continue
          endif

          nin=nin + 1
          input(nin)=line
          istate(nin)=0
        endif
 1500 continue
c
 2000 continue
      write(ierr,*)
      write(ierr,*)'==> Execute'
      return
 2001 write(ierr,*)'Invalid input: EOF; help quits'
      return
c
 3000 write(iout, '(/a)')'                   Rho+ terminated by user'
      stop
 988  write(ierr,*)
     &  'Cannot open '//infil(1:lnblnk(infil))//' to read commands'
      stop 'Rho+ quits'
      end
c______________________________________________________________

      subroutine getone(code, value, nfound)
c$$$$ calls getarr
      implicit double precision (a-h, o-z)
c
c  Extracts a single number from the input store.  That store is
c  a large array in common /store/ which has been filled earlier.
c
c  code    A 4-bye identifying code.  Only lines in the input store
c          beginning with this code are scanned for numbers.
c  value   the real output variable containing the desired number.
c  nfound  is 1 if a number is successfully read in; it is zero
c          the number is absent or unreadable.  nfound = -99 if the
c          code is absent from the input store.
c
      character*4 code
      dimension v(1)
ce
      call getarr(code, v, 1, nfound)
      if (nfound .eq. 1) value=v(1)
      return
      end
c______________________________________________________________

      subroutine getchr(code, cchar, nfound)
c$$$$ calls no other routines
      implicit double precision (a-h, o-z)
c  Extracts a single character variable from the input store.  That 
c  store is a large array in common /store/ which has been filled 
c  earlier.
c
c  code    A 4-byte identifying code.  Only lines in the input store
c          beginning with this code are scanned.
c  char    the character output variable containing the desired string.
c  nfound  is 1 if a string is successfully read in; it is zero if
c          the line is blank after the code.  nfound = -99 if the
c          code is absent from the input store.
c
      parameter (inmx=200)
      character *80 input(inmx),line, cchar*(*),  code*4
      common /store/ input
      common /ndict/ iecho,nin,memory,kread(inmx)
      common /io/ inp,iout,iprint,ierr

c  Inspect the store in reverse order (thus reading latest entry)
ce
      do 1010 lin=nin, 1, -1
        line=input(lin)
c  Check code and if line is fresh or old
        if (code .eq. line(1:4) .and. kread(lin) .le. memory) then
          if (iecho .ge. 1)then
              write(iout,'(a)')' '
              len=min(lnblnk(line),index(line,'%')-1)
              if(len.le.0.0)len=lnblnk(line)
              if(iout.eq.16)
     &           write(iout,'(2a)')'==> '//line(1:len)
              write(ierr,'(a)')' '
              write(ierr,'(2a)')'==> '//line(1:len)
          endif
c  Increment read count
          kread(lin)=kread(lin) + 1
          goto 1020
        endif
 1010 continue
c  Code word not found
      nfound=-99
      return
c
 1020 continue
c  
      n1=index(line, ' ')+1
      n2=index(line, '%')-1
      if (n2 .lt. 0) n2=80
c  Save in cchar everything after 1st blank to just before %
      do 1200 k=n1, n2
        if (line(k:k) .ne. ' ') then
          cchar=line(k:n2)
          nfound=1
          return
        endif
 1200 continue
c
c  Blank field after code word
      nfound=0
      return
      end
c______________________________________________________________

      subroutine getstr(code,string,nfound)
c$$$$ calls no other routines
      implicit double precision (a-h, o-z)
c  Looks for a string argument in a character variable in the
c       input store.  That 
c  store is a large array in common /store/ which has been filled 
c  earlier.
c
c  code    A 4-byte identifying code.  Only lines in the input store
c          beginning with this code are scanned.
c  string  the string argument you are looking for
c  
c  nfound  is index of start of string if it is found
c  nfound  is 0 if the string is not found
c          is -99 if the code is absent from the input store.
c
      parameter (inmx=200)
      character *80 input(inmx),line, string*(*),code*4
      common /store/ input
      common /ndict/ iecho,nin,memory,kread(inmx)
      common /io/ inp,iout,iprint,ierr

c  Inspect the store in reverse order (thus reading latest entry)
ce
      do 1010 lin=nin, 1, -1
        line=input(lin)
c  Check code and if line is fresh or old
        if (code .eq. line(1:4) .and. kread(lin) .le. memory) then

          if (iecho .ge. 1)then
              write(iout,'(a)')' '
              len=min(lnblnk(line),index(line,'%')-1)
              if(len.le.0.0)len=lnblnk(line)
              write(iout,'(2a)')'==> '//line(1:len)
              write(ierr,'(a)')' '
              write(ierr,'(2a)')'==> '//line(1:len)
          endif
c  Increment read count
          kread(lin)=kread(lin) + 1
          goto 1020
        endif
 1010 continue
c  Code word not found
      nfound=-99
      return
c
 1020 continue
c  
      n1=index(line, ' ')+1
      n2=index(line, '%')-1
      if (n2 .lt. 0) n2=80

c  search from just after code to just before % for the
c     the desired string
c
      nfound=ifndst(line(n1:n2),string)
      if(nfound.gt.0)nfound=nfound+n1-1
      return
      end
c______________________________________________________________
 
      integer function ifndst(str1,str2)
      implicit double precision (a-h, o-z)

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
         ifndst=i
         return
 20   continue
 
      ifndst=0
      return
 
      end

c______________________________________________________________________

      subroutine getarr(code, values, nwant, nfound)
c$$$$ calls no other routines
      implicit double precision (a-h, o-z)
c
c  Extracts an array of numbers from the input store.  That store is
c  a large array in common /store/ which has been filled earlier.
c
c  code    A 4-byte identifying code.  Only lines in the input store
c          beginning with this code are scanned for numbers.
c  values  the real output array of values found.
c  nwant   the maximum number of numbers expected .
c  nfound  the number of valid numbers actually found in the input.
c          If the line contains fewer than nwant  values, this is the
c          number returned in  nfound.
c          If there are no numbers after the codeword
c          nfound=0.  Finally, if the code is absent from the store
c          nfound=-99 and the array values is left undisturbed.
c
      parameter (inmx=200)
      character *80 input(inmx),line
      common /store/ input
      common /ndict/ iecho,nin,memory,kread(inmx)
      common /io/ inp,iout,iprint,ierr
      dimension values(*)
      character local*40, code*4, cchar*80
c
c  Read the store in reverse order (Thus last entry is obeyed)
ce
      do 1010 lin=nin, 1, -1
        line=input(lin)
c  Check code and if line is fresh or old
        if (code .eq. line(1:4) .and. kread(lin) .le. memory) then
          if (iecho .ge. 1)then
              write(iout,'(a)')' '
              len=min(lnblnk(line),index(line,'%')-1)
              if(len.le.0.0)len=lnblnk(line)
              write(iout,'(2a)')'==> '//line(1:len)
c             if(iout.eq.16)then
                write(ierr,'(a)')' '
                write(ierr,'(2a)')'==> '//line(1:len)
c             endif
          endif
c  Increment read count; copy tail into cchar, discarding stuff after %
          kread(lin)=kread(lin) + 1
          n1=index(line, ' ')+1
          n2=index(line, '%')-1
          if (n2 .lt. 0) n2=80
          cchar=line(n1:n2)
          kn=n2 - n1 + 1
          goto 1020
        endif
 1010 continue
c  Code word not found
      nfound=-99
      return
c
 1020 continue
      k1=1
c  Up to  nwant  numbers are sought.
      do 1800 nval=1, nwant
        do 1100 k=k1, kn
          if (cchar(k:k) .ne. ' ') goto 1200
 1100   continue
        nfound=nval-1
        return
 1200   do 1300 l=k, kn
          if (cchar(l:l).eq. ',' .or. cchar(l:l).eq.' ') goto 1500
 1300   continue

 1500   if (index('	Ee+-.0123456789', cchar(k:k)) .eq. 0)goto 2000
        if (k-l+1 .gt. 40) goto 2010
        k2=l+1
        local=cchar(k:l-1)
        read(local, '(f40.0)', err=2000) values(nval)
 1800 k1=k2
      nval=1 - nwant
      nfound=1-nval
      return
c
 2000 continue
      nfound=nval-1
      return
 2010 write(iout,'(a)')' ',
     $'>>> Program can not read a number of more than 40 characters:',
     $cchar(k:l-1)
      write(ierr,'(a)')' ',
     $'>>> Program can not read a number of more than 40 characters:',
     $cchar(k:l-1)
      nfound=1-nval
      return
      end
c_______________________________________________________________________

      subroutine report()
c$$$$ calls mult, forwrd, span
c  Prints a solution summary, giving transition points for vector
      implicit double precision (a-h, o-z)
      parameter (mxf=101,mxxi=450, mxrw=2*mxf)
      common // a(mxrw*mxxi),d(mxrw),x(mxxi),ma,na
      common /work/ w(mxxi),zz(mxrw),act(mxrw,mxrw+2),istt(mxxi+1)
      common /observ/ freq(mxf),rho(mxf),drh(mxf),phi(mxf),dph(mxf)
     $  ,erh(mxf),eph(mxf),dwtrh(mxf),dwtph(mxf),sgr(mxf),nfrq
      common /shared/ pi,twopi,amu0,radian
      common /discrt/ xi(mxxi),nxi
      common /final/ respon(2,mxf),bfresp(2,mxf),bound(2,2*mxf)
      common /prolix/ loquor
      logical logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs(4)
      character*3 cecho
      common /flags/ logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs,cecho
      common /io/ inp,iout,iprint,ierr
      character*40 part1,part2
c
      if(loquor.gt.0)write(iout,'(/a,1p,g12.5)')' ln A: ',x(1)
      if (loquor .ge. 2) then
        write(iout,'(a)')' ','  Transitions in solution vector',
     $      '          lambda   mu(n-1) mu(n) mu(n+1)'
        j1=1
        do 1100 j=3, nxi
          if (x(j) .ne. x(j-1) .and. j.gt. j1) then
            j1=min(nxi, j+1)
            write(iout,'(i7,1p,g12.4,0p,3f6.2)')
     $          j,xi(j-1),x(j-1),x(j),x(j1)
          endif
 1100   continue
      endif
c
c
c  Print residuals in internal format - debugging only
      if (loquor .ge. 4) then
c  Find predictions of current model, setting them in w
        call mult (ma, nxi, 1, a,x, w)
        i=0
        write(9,'(a)')' ',
     $  '           Table of misfits to model predictions',
     $  '      freq         Re dat    Re pred    diff     '
     $                   //'Im dat    Im pred    diff'
        do 1200 j=1, nfrq
          part1='       --         --        -'
          if (drh(j) .eq. 0 .and. dph(j) .eq. 0) goto 1200
          if (drh(j) .gt. 0) then
            i=i+1
            write(part1, '(2f11.1,f8.1)') d(i),w(i),w(i)-d(i)
          endif
          part2='       --         --        -'
          if (dph(j) .gt. 0) then
            i=i+1
            write(part2, '(2f11.1,f8.1)') d(i),w(i),w(i)-d(i)
          endif
          write(9,'(i3,1pg12.4,2a)')j,freq(j),part1,part2
 1200   continue
      endif
c
c  Find model predictions at all freqs - results in /final/ .. respon
      call forwrd(nfrq, freq, respon)

      if (loquor .ge. 0) then
         if(logdat)then
           write(iout,'(a)')' ',
     $'            Table of misfits between data and model predictions',
     $' ',
     $'                   - Log10 Apparent Resistivity -   '
     $                  //'---------- Phase ----------',
     $'    freq           data    predicted  %dR/R  norm   '
     $                  //'data  predicted   dP   norm',
     $'                                                    '
     $                  //'               (scaled)     '
         else
           write(iout,'(a)')' ',
     $'            Table of misfits between data and model predictions',
     $' ',
     $'                ----- Apparent Resistivity ------   '
     $                  //'---------- Phase ----------',
     $'    freq        data    predicted     %dR/R  norm   '
     $                  //'data  predicted   dP   norm',
     $'                                                    '
     $                  //'               (scaled)    '


         endif
      endif

      do 1400 j=1, nfrq
        part1='     --        --           -     -     '
        if(logdat)then
          write(part1(12:22),'(f11.4)') log10(respon(1,j))
        else
          write(part1(12:22),'(1p,g11.4)') respon(1,j)
        endif

        if (drh(j) .gt. 0) then
          pred=respon(1,j)
          if(logdat)then
            write(part1,'(2f11.4,f8.1,f6.1)') log10(rho(j)),
     &                                     log10(pred),
     &                             230.29*(log10(pred)-log10(rho(j))),
     &                    rho(j)*abs(log(pred)-log(rho(j)))/drh(j)
          else
            write(part1,'(1p,2g11.4,0p,f8.1,f6.1)') rho(j),
     &                                     pred,
     &                                   100.*(pred-rho(j))/rho(j),
     &                    rho(j)*abs(log(pred)-log(rho(j)))/drh(j)
          endif
        endif

        part2='    --      --       -      -'

        write(part2(9:16),'(f8.2)') respon(2,j)

        if (dph(j) .gt. 0) then
          pred=respon(2,j)
          write(part2, '(2f8.2,f8.1,f6.1)') phi(j),
     &                                      pred,
     &                                      (pred-phi(j))*3.4907,
     &                                      abs(pred-phi(j))/dph(j)
        endif

        if (loquor .ge. 0) write(iout,'(i3,1pg10.3,a36,a30)')
     $  j,freq(j),part1,part2
c
 1400 continue

        if (loquor .ge. 0) then
         write(iout,*)
         write(iout,*)
     &     ' - norm: missfit divided by data error'
     &     //' (as modified by downweights)'
         if(logdat)then
           write(iout,*)
     &       ' - %dR/R: (log10(rho) misfit) * 230.29 ='
     &       //' (resistivity missfit)/resistivity in %'
         else
           write(iout,*)
     &       ' - %dR/R: (rho missfit)/rho in %'
         endif
         write(iout,*)
     &     ' - dP(scaled): (phase missfit) * 3.491 to give it'
     &     //' same scale as %dR/R'
         write(iout,*)
        endif

      return
      end
c_______________________________________________________________________

      subroutine printa(title, mdim, m ,n, a)
c$$$$ calls no other routines
c  Print a vector or an array: handy debugging routine

      implicit double precision (a-h, o-z)
      common /io/ inp,iout,iprint,ierr
      dimension a(mdim,*)
      character*(*) title
c
      write(iout,*)title
      if (n.eq. 1) then
        write(iout,'(1p6g12.4)')(a(j,1),j=1, m)
      else
        do 1100 i=1, m
          if (m.gt.1) write(iout,*)'Row',i
          write(iout,'(1p6g12.4)')(a(i,j),j=1, n)
 1100   continue
      endif
      return
      end
c______________________________________________________________________

      subroutine plot(pltall)
c$$$$ calls getchr

c  calls a user-defined plotting routine

      implicit double precision (a-h, o-z)

      logical pltall

c generates a plot of only the data and the rho+ responses if
c     pltall is .false. Otherwise generates full plot including
c     any bounds.

      parameter (mxf=101,mxxi=450, mxrw=2*mxf)
      parameter (inmx=200)
      common // a(mxrw*mxxi),d(mxrw),x(mxxi),ma,na
      common /work/ w(mxxi),zz(mxrw),act(mxrw,mxrw+2),istt(mxxi+1)
      common /observ/ freq(mxf),rho(mxf),drh(mxf),phi(mxf),dph(mxf)
     $  ,erh(mxf),eph(mxf),dwtrh(mxf),dwtph(mxf),sgr(mxf),nfrq
      common /shared/ pi,twopi,amu0,radian
      common /final/ respon(2,mxf),bfresp(2,mxf),bound(2,2*mxf)
      common /chore/ freqs(2*mxf),ntask
      common /discrt/ xi(mxxi),nxi
      common /prolix/ loquor
      common /optiml/ wtmax,fitmin,prob
      common /io/ inp,iout,iprint,ierr
      logical logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs(4)
      character*3 cecho
      common /flags/ logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs,cecho
      common /ndict/ iecho,nin,memory,istate(inmx)

      character*75 pfile, root*70
      character*70 pltype

c  Unit number for plot output:
      iplt=10

c  Determine if plot requested and what kind:
      call getchr('plot', pltype, iplot)
      if (iplot .lt. 0)then
          return
      else
c  If plot requested, get root name for files if it exists
        if(cecho.eq.'on')iecho=-1
        call getchr('root', root, isroot)
        if(cecho.eq.'on')iecho=1
        if (isroot .gt. 0)then
          length=index(root,' ')-1
          pfile=root(1:length)
          pfile(length+1:length+4)='.plt'
        else
          pfile='rplus.plt'
        endif
        open (unit=iplt, file=pfile, err=999)

c  If plot requested, call routine to make plotfile

        if(pltype.eq.'plotxy')then
           call plt_uw(pltall,iplt,.true.,pfile)
        elseif(pltype.eq.'xyplot')then
           call plt_uw(pltall,iplt,.false.,pfile)

CCC@
c     implement additional choices here:
         elseif(pltype.eq.'myplot')then
c         user may want to edit:
            call myplot(pltall,iplt)
        else

c  Default if plot type illegal or unspecified (user should edit)
           pltype='xyplot'
           call plt_uw(pltall,iplt,.false.,pfile)
        endif

        close(unit=iplt)

        write(iout,'(a)')
     &    ' ','Plotfile type is: '//pltype(1:lnblnk(pltype))
        write(iout,'(a)')
     &    ' ','Plotfile written to '//pfile(1:lnblnk(pfile))
        write(iout,*)' '

        return
      endif

999   write(ierr,*)'No plot file. Cannot open '//pfile(1:lnblnk(pfile))
      write(iout,*)'No plot file. Cannot open '//pfile(1:lnblnk(pfile))
      return

      end

c______________________________________________________________________

      subroutine myplot(pltall,iplt)
 
CCC@
c  user-defined plotting routine
 
      implicit double precision (a-h, o-z)

c  input arguments:
      integer iplt
c        unit number of plot output file
c      integer isroot
c        >0 root name for plotfile given by user
c     character*(*) pfile
c         name of plot file
      logical pltall
c         if .false. plot only data and responses
c         if .true.  plot bounds as well as data and responses

c     iplt is opened in calling routine, thus isroot and pfile
c       are provided for information purposes only and can be
c       used for plot labels, etc.


c  commons to import numbers to plot:

      parameter (mxf=101,mxxi=450, mxrw=2*mxf)
      parameter (inmx=200)
      common // a(mxrw*mxxi),d(mxrw),x(mxxi),ma,na
      common /work/ w(mxxi),zz(mxrw),act(mxrw,mxrw+2),istt(mxxi+1)
      common /observ/ freq(mxf),rho(mxf),drh(mxf),phi(mxf),dph(mxf)
     $  ,erh(mxf),eph(mxf),dwtrh(mxf),dwtph(mxf),sgr(mxf),nfrq
      common /shared/ pi,twopi,amu0,radian
      common /final/ respon(2,mxf),bfresp(2,mxf),bound(2,2*mxf)
      common /chore/ freqs(2*mxf),ntask
      common /discrt/ xi(mxxi),nxi
      common /prolix/ loquor
      common /io/ inp,iout,iprint,ierr
      common /optiml/ wtmax,fitmin,prob
      logical logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs(4)
      character*3 cecho
      common /flags/ logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs,cecho
      common /ndict/ iecho,nin,memory,istate(inmx)
c 
c Locations and definitions of numbers likely to be used in a plot:
c
c  References to data starts at index 2 because all data, errors
c    and weights were shifted by 1 in subroutine range
c
c        freq(j)  j=2,nfrq: data frequencies
c        erh(j)   j=2,nfrq: original errorbar for rho
c        eph(j)   j=2,nfrq: original errorbar for phi
c        dwtrh(j) j=2,nfrq: downwight factor for rho (=0 if excluded)
c        dwtph(j) j=2,nfrq: downwight factor for phi (=0 if excluded)
c        drh(j)   j=2,nfrq: weighted errorbar that is used for inversion
c        dph(j)   j=2,nfrq: weighted errorbar that is used for inversion
c
c  References to final responses and bounds are NOT shifted like data
c
c        bfresp(1,j) j=1,nfrq: Rho-a of best-fitting model in D+
c        bfresp(2,j) j=1,nfrq: Phase of best-fitting model in D+
c
c        prob*100+.5 is confidence level for bounds
c
c        freqs(j) j=1,ntask: frequencies at which bounds are computed
c         NOTE that freqs<0 correspond to phase bounds
c         NOTE that freqs>0 correspond to rho-a bounds
c        bound(1,j) j=1,ntask AND freqs(j)<0: phase bounds
c        bound(2,j) j=1,ntask AND freqs(j)>0: rho-a bounds
c

c     put plotting code here to write plot commands for data
c        and calculated responses to unit iplt

      write(iplt,*)'myplot not implemented'
      write(iplt,*)'   No data or response plot commands written'

      if(pltall)then
c        put plotting code here to write plot commands for
c           bounds to unit iplt
         write(iplt,*)'   No bound plot commands written'
      endif

      return
      end
c______________________________________________________________________

      subroutine plt_uw(pltall,iplt,ucsd,pfile)
c$$$$ calls span
      implicit double precision (a-h, o-z)

      integer iplt
      logical ucsd, pltall
      character*(*) pfile

c  If pltall=.true. Plots all numbers (data, responses, bounds)
c  If pltall=.false. Plots only data and responses

c  If ucsd=.true., writes an input file for Parker's plotxy program
c     summarizing results.
c  If ucsd=.false. file is modified to work with xyplot, which is
c   the U.W./Chave/Schultz version of plotxy.
c
c  The file will generate graphs of:
c
c  Original data:
c     Included data are shown as solid small circles.
c     Excluded data are shown as open small circles.
c     Error bars are modified by the weights requested.
c     Rho-a's =< 0 are never plotted.
c
c  Responses of the best-fitting model at all frequencies in
c     the data file are shown as open octagons
c
c  Bounds requested are shown as triangles whose apexes point
c    up for lower bounds and down for upper bounds. These
c    points are connected with dashed lines.
c  NOTE: these dashed lines are for convenience only and have NO
c    OTHER SIGNIFICANCE

c  Information comes through /observ/, /final/, /optiml/,/chore/
      parameter (mxf=101,mxxi=450, mxrw=2*mxf)
      parameter (inmx=200)
      common // a(mxrw*mxxi),d(mxrw),x(mxxi),ma,na
      common /work/ w(mxxi),zz(mxrw),act(mxrw,mxrw+2),istt(mxxi+1)
      common /observ/ freq(mxf),rho(mxf),drh(mxf),phi(mxf),dph(mxf)
     $  ,erh(mxf),eph(mxf),dwtrh(mxf),dwtph(mxf),sgr(mxf),nfrq
      common /shared/ pi,twopi,amu0,radian
      common /final/ respon(2,mxf),bfresp(2,mxf),bound(2,2*mxf)
      common /chore/ freqs(2*mxf),ntask
      common /discrt/ xi(mxxi),nxi
      common /prolix/ loquor
      common /io/ inp,iout,iprint,ierr
      common /optiml/ wtmax,fitmin,prob
      logical logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs(4)
      character*3 cecho
      common /flags/ logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs,cecho
      common /ndict/ iecho,nin,memory,istate(inmx)

c Setting fancy=.false. reduces plot complexity
      logical fancy

      if(cecho.eq.'on')iecho=-1
      call getstr('plot','fan',nfound)
      if(cecho.eq.'on')iecho=1
      if(nfound.gt.0)then
         fancy=.true.
      else
         fancy=.false.
      endif
c
c  Convert frequency to period for plotting 
c  (IAGA standard requires long period on right of plot and
c   this is the easiest way to obey the standard.)

c  References to data starts at index 2 because all data, errors
c    and weights were shifted by 1 in subroutine range
      do 10 j=2,nfrq
10      freq(j)=1./freq(j)
      do 12 j=1,ntask
12      freqs(j)=1./freqs(j)

c  Set up initial state of plotx
      if (ucsd) then
        write(iplt,'(a)') 'frame on'
      else 
        write(iplt,'(a)') 'frame 1','stack 1'
      endif
      write(iplt,'(a)') 'file *','logxy 1',
     $'xlim 6 ','ylim 4',
     $'xlabel Period (s)','ylab Phase (degrees)'
c
c  Count number of phase and rhoa points required
c   (also count number of rhoa<=0 points to be excluded entirely)
      np=0
      nr=0
      nrn=0
c  Starts at 2 because 1st row is a bounding constraint
c  (All data were shifted by one in subroutine range)
      do 1100 j=2, nfrq
        if (dwtph(j) .gt. 0.0) np=np+1
        if (dwtrh(j) .gt. 0.0 .and. erh(j).ge. 0.0) nr=nr+1
        if (erh(j) .lt. 0.0) nrn=nrn+1
 1100 continue
c
c  First plot - Phi information

c  Included orginal data with uncertainties as used in calculation.
      if (np .gt. 0)then
          if (.not. ucsd) write(iplt,'(a)') 'smooth -1'
          write(iplt,'(a)')'symbol 19','mode 3'
          write(iplt,'(a,i4)')'read',np
      endif
      do 1200 j=2, nfrq
        if (dwtph(j) .gt. 0.0) write(iplt,'(1p,3g12.4)')
     $  freq(j), phi(j), dph(j)
 1200 continue

c  Excluded orginal data with uncertainties
c  Excluded data (as modified by gtband) with uncertainties
      if(nfrq-np-1.gt.0)then
         if (.not. ucsd) write(iplt,'(a)') 'smooth -1'
         write(iplt,'(a)')'symbol 16 0.25','mode 3'
         write(iplt,'(a,i4)')'read',nfrq-1-np
         do 1210 j=2, nfrq
           if (dwtph(j) .eq. 0.0) write(iplt,'(1p,3g12.4)')
     $     freq(j), phi(j), eph(j)
 1210    continue
      endif

c  Calculated bound values
      if(pltall)then
        np=0
        do 1300 nt=1, ntask
          if (freqs(nt) .lt. 0.0) np=np + 1
 1300   continue
c
        if (np .gt. 0) then
          if(fancy)then
              write(iplt,'(a)') 'smooth -1'
              write(iplt,'(a)')'symbol 2 .09'
          else
              if (.not. ucsd) write(iplt,'(a)') 'smooth 0'
              write(iplt,'(a)')'dash .03 .03'
          endif
          write(iplt,'(a)')'mode 2'
          write(iplt,'(a,i4)')'read',np

          do 1400 j=1, ntask
           if (freqs(j).lt.0.0)
     &        write(iplt,'(2g12.4)')-freqs(j),bound(1,j)
 1400     continue
          if(fancy)then
              write(iplt,'(a)') 'smooth -1'
              write(iplt,'(a)')'symbol 21'
          else
             if (.not. ucsd) write(iplt,'(a)') 'smooth 0'
             write(iplt,'(a)')'dash .03 .03'
          endif
          write(iplt,'(a)')'mode 2'
          write(iplt,'(a,i4)')'read',np

          do 1500 j=1, ntask
            if (freqs(j).lt.0.0)
     &         write(iplt,'(2g12.4)')-freqs(j),bound(2,j)
 1500     continue
        endif
      endif

c  Response of best-fitting model (response NOT shifted like data)
      if(.not. ucsd)write(iplt,'(a)')'smooth 0'
      write(iplt,'(a)')'mode 2'
      write(iplt,'(a,i4/(2g12.4))')'read',nfrq-1,
     $(freq(j),bfresp(2,j-1),j=2, nfrq)
      if (ucsd) then
        write(iplt,'(a)') 'plot 1.5 1.0', 'stack'
      else
        write(iplt,'(a)') 'offset 1.5 1.0','plot'
      endif
c
c  Second plot - Rho-a information

      write(iplt,'(a)') 'xlab','ylab Log10 \\r\\'// '\\sub{a} (ohm m)'
c     write(iplt,'(a)') 'xlab','ylab \\r\\'// '\\sub{a} (ohm m)'

c     if (isroot .lt. 0)then
        if(pltall)then
          write(iplt,'(a9,1pg11.4,a4,i3,a2,a10)')'title Xsq',
     $     fitmin,' Bnd',
     $     int(prob*100+.5),'% ',pfile(1:lnblnk(pfile)-4)
        else
          write(iplt,'(a16,1pg11.4,a,a12)')'title Rho+ Chisq',
     $     fitmin,' ',pfile(1:lnblnk(pfile)-4)
        endif
c     else
c       if(pltall)then
c         write(iplt,'(a9,1pg11.4,a4,i3,a2,a10)')'title Xsq',
c    $     fitmin,' Bnd',
c    $     int(prob*100+.5),'% ',pfile(1:lnblnk(pfile)-4)
c       else
c         write(iplt,'(a16,1pg11.4,a,a12)')'title Rho+ Chisq',
c    $     fitmin,' ',pfile(1:lnblnk(pfile)-4)
c       endif
c     endif
c
c     call span(nfrq, rho, drh, r1, r2)
c
c  Included original data with uncertainties as used by calculation
      write(iplt,'(a)')'ylim 4'
      if (.not. ucsd) write(iplt,'(a)') 'smooth -1'
      write(iplt,'(a)')'symbol 19','logxy 1'
c     write(iplt,'(a)')'symbol 19','logxy 3'
      write(iplt,'(a/a,i5)') 'mode 3','read',nr
      do 2200 j=2, nfrq
        if (dwtrh(j) .gt. 0.0) write(iplt,'(1p,3g12.4)')
     $  freq(j), log10(rho(j)), .4343*drh(j)/rho(j)
c    $  freq(j), rho(j), drh(j)
 2200 continue
      
c  Excluded original data with uncertainties (except rho-a =< 0)
c  Excluded data (as modified by gtband) with uncertainties
      if(nfrq-nr-nrn-1 .ne. 0)then
         if (.not. ucsd) write(iplt,'(a)') 'smooth -1'
         write(iplt,'(a)')'symbol 16 .25','logxy 1'
c        write(iplt,'(a)')'symbol 16 .25','logxy 3'
         write(iplt,'(a/a,i5)') 'mode 3','read',nfrq-1-nr-nrn
         do 2210 j=2, nfrq
           if (dwtrh(j) .eq. 0.0 .and. erh(j).ge.0.0)
     $       write(iplt,'(1p,3g12.4)')
     $         freq(j), log10(rho(j)), .4343*erh(j)/rho(j)
c    $         freq(j), rho(j), erh(j)
 2210    continue
      endif

c  Calculated bound values
      if(pltall)then
        nr=ntask-np
        if (nr .gt. 0) then
          if(fancy)then
             write(iplt,'(a)') 'smooth -1'
             write(iplt,'(a)')'symbol 2 .09'
          else
             if (.not. ucsd) write(iplt,'(a)') 'smooth 0'
             write(iplt,'(a)')'dash .03 .03'
          endif
          write(iplt,'(a)')'mode 2'
          write(iplt,'(a,i4)')'read',nr

          do 2400 j=1, ntask
           if (freqs(j).gt.0.0)
     $        write(iplt,'(2g12.4)')freqs(j),log10(bound(1,j))
c    $        write(iplt,'(2g12.4)')freqs(j),bound(1,j)
 2400     continue

          if(fancy)then
             write(iplt,'(a)') 'smooth -1'
             write(iplt,'(a)')'symbol 21'
          else
             if (.not. ucsd) write(iplt,'(a)') 'smooth 0'
             write(iplt,'(a)')'dash .03 .03'
          endif
          write(iplt,'(a)')'mode 2'
          write(iplt,'(a,i4)')'read',nr

          do 2500 j=1, ntask
           if (freqs(j).gt.0.0)write(iplt,'(2g12.4)')
     $         freqs(j),log10(bound(2,j))
c    $         freqs(j),bound(2,j)
 2500     continue
        endif
      endif
c
c  Response of best-fitting model (response NOT shifted like data)
      write(iplt,'(a)')'mode 2'
      write(iplt,'(a)')'smooth 0'
      write(iplt,'(a,i4/(2g12.4))')'read',nfrq-1,
     $(freq(j),log10(bfresp(1,j-1)),j=2, nfrq)
c    $(freq(j),bfresp(1,j-1),j=2, nfrq)
c
      if (.not. ucsd) write(iplt,'(a)') 'offset 0 4.5'
      write(iplt,'(a)')'plot','stop'
      return
      end

c_______________________________________________________________________

      subroutine gtband()
c$$$$ calls getchr, getarr, ckband

      implicit double precision (a-h, o-z)

c this subroutine communicates with the outside through

      parameter (inmx=200)
      common /ndict/ iecho,nin,memory,istate(inmx)
      common /bopram/ iband,nbands,nstart(4),nend(4),gain(4)
     
      logical logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs(4)
      character*3 cecho
      common /flags/ logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs,cecho
      common /io/ inp,iout,iprint,ierr
 
c nstart and nend are arrays of the data indices of the band edges.
c nstart must be < nend for each band. The bands are incorporated
c in the calculation in the order read in.
c
c See function ckband for other limitations wrt to chosen bands
 
c A maximum of 4 non-overlapping bands are allowed and the first
c   one read is assumed to have no offset.
c
c Subroutine fndofs holds the first band read in fixed and computes
c   the band offsets in the order read in.
c  
c The comand key word to invoke band offsets is "band" followed
c  by the band edge indices (bottom then top) and an optional
c  option number (0 to 5). 

c The option controls setting of the following flags:

c  phjudg:  false - use chisq phi to judge offset
c           true  - use chisq to judge offset
c  usephi:  true  - phase data used in computing band offsets
c           false - phase data not used in computing band offsets
c  phfit:   true  - use phase data in computing responses and misfit
c                   after rho bands are offset

c The flag "rhoflg": can be changed in this routine at compile time.
c   Experience has shown that it makes little difference what its
c   value is as long as usephi=true

c  rhoflg:  true  - Rho data, in bands whose offset has already
c                   been computed, are included when the next band
c                   offset is computed.
c           false - Only the fixed band and the band, whose offset is
c                   being computed are included during each offset
c                   calculation. (Note this option should be used
c                   with great care, especially if usephi=false).

c  ioptn: missing : offst=F ioptn = 0 (DEFAULT)
c  ioptn: illegal : offst=F (offsets will not be computed, BUT may
c                    be applied by user using the gain command
c
c  ioptn=-1 to -5 : compute band offsets:
c
c        -1 : phjudg=F; usephi=T; phfit=T  Use all data all the time.
c        -2 : phjudg=F; usephi=F; phfit=F  Use only rho data for gains and
c                                           any subsequent tasks.
c                                           (this is equivalent to setting
c                                            all phi dwnwts to 0.0)
c        -3 : phjudg=F; usephi=F; phfit=T  Use only rho data for gains, but
c                                           include phase in final fit
c                                           and any subsequent tasks.
c Following options are not recommended for general use. -5 is particularly
c    dangerous and should be used with great caution.
c        -4 : phjudg=T; usephi=T; phfit=T  Use all data, but judge
c                                            gain by phase fit only.
c        -5 : phjudg=T; usephi=F; phfit=T  Use only rho data to determine
c                                            best-fitting model at each
c                                            gain iteration, but judge
c                                            gain by the phase fit of
c                                            each such model.
c                                          Phase included in final fit
c                                           and any subsequent tasks.

      dimension wrk(8)
      character*40 chnull

      call getchr('band',chnull,need)
      if(need.gt.0)then

        if(cecho.eq.'on')iecho=-1
        call getarr('band',wrk,9,nbands)
        if(cecho.eq.'on')iecho=1

        if(mod(nbands,2).gt.0 .and. nbands.gt.1)then
c  there is an option flag
           if(wrk(nbands).le.-1. .and. wrk(nbands).ge.-5.)then
             ioptn= wrk(nbands)
           else
c  option flag value out of range (-5,-1), default is
             ioptn=0
           endif
        else
c  there is no legal option flag, default is
           ioptn= 0
        endif

        rhoflg=.true.

        if(ioptn.eq.-1)then
          phjudg=.false.
          usephi=.true.
          phfit=.true.
        elseif(ioptn.eq.-2)then
          phjudg=.false.
          usephi=.false.
          phfit=.false.
        elseif(ioptn.eq.-3)then
          phjudg=.false.
          usephi=.false.
          phfit=.true.
        elseif(ioptn.eq.-4)then
          phjudg=.true.
          usephi=.true.
          phfit=.true.
        elseif(ioptn.eq.-5)then
          phjudg=.true.
          usephi=.false.
          phfit=.true.
        else
          phjudg=.false.
          usephi=.true.
          phfit=.true.
        endif

        nbands=nbands/2
        if(ioptn.ge.0 .or. nbands.lt.2)then
           offst=.false.
        else
           offst=.true.
        endif

        if(nbands.gt.0)then
          do 30 i=1,nbands
            nstart(i)=wrk(2*i-1)
            nend(i)=wrk(2*i)
 30       continue
        else
          return
        endif
      else
        nbands=0
        offst=.false.
        phfit=.true.
        return
      endif
 
c test the requested bands to see if they are compatible with the
c  input data.
 
      call ckband(ioptn)
      if(.not. flg2) nbands=0

c  see if there are any fixed band offsets to be applied prior
c   to any computations

      call getchr('gain',chnull,need)
      if(need.gt.0)then
        if(cecho.eq.'on')iecho=-1
        call getarr('gain',wrk,4,n)
        if(cecho.eq.'on')iecho=1

        if(n.ne.nbands)then
          write(ierr,*)' '
          write(ierr,*)'>>> Warning: # of gains not equal to # of bands'
          write(ierr,*)'    No gains applied prior to computation'
          write(iout,*)
          write(iout,*)'>>> Warning: # of gains not equal to # of bands'
          write(iout,*)'    No gains applied prior to computation'
          fixgn=.false.
        else
          fixgn=.true.
          do 100 i=1,nbands
             if(wrk(i).ne.1.0)then
                if(logdat)then
                  gain(i)=wrk(i)
                else
                  gain(i)=log10(wrk(i))
                endif
                incofs(i)=.false.
c               incofs(i)=.true.
             else
                gain(i)=0.0
                incofs(i)=.true.
             endif
 100      continue
          endif
      else
        do 101 i=1,nbands
          gain(i)=0.0
          incofs(i)=.true.
          fixgn=.false.
 101    continue
      endif
 
      return
      end
 
c_______________________________________________________________________
 
      subroutine ckband(ioptn)
      implicit double precision (a-h, o-z)
c$$$$$ calls indxnt
 
c returns flg2 (in common /flags/ .true. if band offsets are compatible
 
c Compatibility requires:
c   Band start index must be <= band end index
c   Non-overlapping bands
c   At least one included resistivity in each band unless offsets
c     are not being calculated in which case there must be at least
c     one included rho in some band.
c   At least one included phase datum, unless ioptn = -3
c   Bands must not extend beyond data

c Addtionally, this routine changes the included data in the following ways:
c   If ioptn>=0 (offsets not calculated) then ALL data outside of
c       the specified bands are excluded.
c   If ioptn<0 (offsets calculated) rho data outside of
c       the specified bands are excluded. Exclusion of phase data outside
c       bands requirs the prunephase option or
c       setting their downweights to zero.
 
      parameter (mxf=101,mxxi=450, mxrw=2*mxf)
      parameter (inmx=200)

      common /observ/ freq(mxf),rho(mxf),drh(mxf),phi(mxf),dph(mxf)
     $  ,erh(mxf),eph(mxf),dwtrh(mxf),dwtph(mxf),sgr(mxf),nfrq
      common /bopram/ iband,nbands,nstart(4),nend(4),gain(4)
      logical logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs(4)
      common /ndict/ iecho,nin,memory,istate(inmx)
      character*3 cecho
      common /flags/ logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs,cecho
      common /io/ inp,iout,iprint,ierr

      external Vsum

      integer iordst(4), iordnd(4)
      logical prune, rhprn, phprn

c     write(0,*)'nbands ',nbands
c     write(0,*)'nstart ',(nstart(i),i=1,nbands)
c     write(0,*)'nend ',(nend(i),i=1,nbands)

c Make sure there is at least one frequency in each band

      do 10 i=1,nbands
        if(nend(i).lt.nstart(i))then
          write(ierr,*)'Band ',i,': end index < start index'
          stop '    Rho+ quits'
        endif
 10   continue

c create index matrix pointing to nstart and nend in increasing order

      call indxnt(nstart,nbands,iordst)
c     write(0,*)'iordst ',(iordst(i),i=1,nbands)
      call indxnt(nend,nbands,iordnd)
c     write(0,*)'iordnd ',(iordnd(i),i=1,nbands)

c band ends must not extend outside of data

      do 15 i=1,nbands
        if(nstart(i).lt.1)then
           write(ierr,*)'Band ',i,': extends below data'
           stop '    Rho+ quits'
        endif
        if(nend(i).gt.nfrq)then
           write(ierr,*)'Band ',i,': extends above data'
           stop '    Rho+ quits'
        endif
 15   continue

c bands must not overlap

      if(nbands.gt.1)then
        do 20 i=1,nbands-1
          if(nstart(iordst(i+1)) .le. nend(iordnd(i)))then
            write(ierr,*)'Band ',iordst(i),' overlaps band ',iordnd(i+1)
            stop '    Rho+ quits'
          endif
 20     continue
      endif

c iordst and iordnd must agree

      do 30 i=1,nbands
        if(iordst(i).ne.iordnd(i))then
          write(ierr,*)
     &      'Order of band start and end indices do not agree'
          write(ierr,*)'  start index order: ',(iordst(j),j=1,nbands)
          write(ierr,*)'  end index order  : ',(iordnd(j),j=1,nbands)
          stop '    Rho+ quits'
        endif
 30   continue

c Prune data outside of bands.

c *** NOTE if offst=.true., phase will be pruned only if prune=.true.

      if(cecho.eq.'on')iecho=-1
      call getstr('band','p',nfound)
      if(cecho.eq.'on')iecho=1
      if(nfound.gt.0)then
        prune=.true.
      else
        prune=.false.
      endif

c     write(0,*)'prune ',prune
c     write(0,*)'offst ',offst

c Prune data below lowest band

      if(nstart(iordst(1)).gt.1)then
        do 40 i=nstart(iordst(1))-1,1,-1
          rhprn=.true.
          drh(i)=0.
          dwtrh(i)=0.
          if(.not. offst .or. prune)then
             phprn=.true.
             dph(i)=0.
             dwtph(i)=0.
          endif
 40     continue
      endif
C
c Prune data between bands

      do 55 j=2,nbands
        if(nend(iordnd(j-1))+1.lt.nstart(iordst(j)))then
          do 50 i=nend(iordnd(j-1))+1,nstart(iordst(j))-1
            rhprn=.true.
            drh(i)=0.
            dwtrh(i)=0.
            if(.not. offst .or. prune)then
               phprn=.true.
               dph(i)=0.
               dwtph(i)=0.
            endif
 50       continue
        endif
 55   continue

c Prune data above highest band

      if(nend(iordnd(nbands)).lt.nfrq)then
        do 60 i=nend(iordnd(nbands))+1,nfrq
          rhprn=.true.
          drh(i)=0.
          dwtrh(i)=0.
          if(.not. offst .or. prune)then
             phprn=.true.
             dph(i)=0.
             dwtph(i)=0.
          endif
 60     continue
      endif

      if(rhprn .or. phprn) then
           write(iout,*)
           write(ierr,*)
      endif
      if(rhprn)then
         write(iout,*)'Rho data excluded outside bands'
         write(ierr,*)'Rho data excluded outside bands'
         if(phprn)then
           write(iout,*)'Phase data excluded outside bands'
           write(ierr,*)'Phase data excluded outside bands'
         else
           write(iout,*)'Phase data not pruned outside bands'
           write(ierr,*)'Phase data not pruned outside bands'
         endif
      endif

c If phase are used in offset calculation, make sure there is at least one
c   included phase datum
      if(ioptn .eq. -2 .or. ioptn .eq. -5)then
        sum=Vsum(dph,nfrq)
        if(sum.eq.0.)then
          write(ierr,*)'Offset option ',ioptn,' requires phase data.'
     &     //' No phase included.'
          stop '     Rho+ Quits'
        endif
      endif
 
c If offst=.true. make sure there is at least one rho in each band,
c    otherwise, make sure there is still one included rho

      if(offst)then
       do 120 i=1,nbands
         n=nend(iordnd(i))-nstart(iordst(i))+1
         sum=Vsum( drh(nstart(iordst(i))),n)
         if(sum .eq. 0.)then
           write(ierr,*)'Band ',i,' has no included rho'
           stop
     &     'Must be a rho in each band to calculate offsets. Rho+ quits'
         endif
 120   continue
      else
        sum=Vsum(drh,nfrq)
        if(sum.eq.0.)then
           stop 'No rho included in any band. Rho+ quits'
        endif
      endif

c If you get here, bands are compatible

      flg2=.true.
      return
      end
c_______________________________________________________________________

      function Vsum(x,n)
      implicit double precision (a-h, o-z)
c     Vsum()= sum x(i)
      integer i,n
      dimension x(n)
 
      Vsum=0.0
      do 100 i=1,n
          Vsum= Vsum + x(i)
 100  continue
      return
      end

c_______________________________________________________________________

      subroutine indxnt(ix,n,indx)
      integer n, indx(n)
      integer ix(n)
c      Creates an index, indx(), of the array ix(), referencing ix()
c      in ascending order.
c      Uses a heapsort algorthim to order the references, leaving ix()
c      untouched.
c      integer version of indexx from Torquil Smith's mvlib

      integer i, indxt, ir, j, l
      integer ixt

      do 11 j=1,n
          indx(j)=j
11    continue

      if(n.eq.1) return

      l=n/2+1
      ir=n
10    continue
c         indent
          if(l.gt.1)then
              l=l-1
              indxt=indx(l)
              ixt=ix(indxt)
          else
              indxt=indx(ir)
              ixt=ix(indxt)
              indx(ir)=indx(1)
              ir=ir-1
              if(ir.eq.1)then
                  indx(1)=indxt
                  return
              endif
          endif

          i=l
          j=l+l
20        if(j.le.ir)then
c                 indent
                  if(j.lt.ir)then
                      if( ix(indx(j)).lt.ix(indx(j+1)) )j=j+1
                  endif
                  if( ixt.lt.ix(indx(j)) )then
                      indx(i)=indx(j)
                      i=j
                      j=j+j
                  else
                      j=ir+1
                  endif
c                 unindent
              go to 20
          endif

          indx(i)=indxt
c         unindent
      go to 10
      end

c_______________________________________________________________________

      subroutine bvls(key, m, n, a, b, bl, bu, x, w, act, zz, istate,
     +  loopA)
      implicit double precision (a-h, o-z)
c
c$$$$ calls qr
c--------------------Bounded Variable Least Squares---------------------
c
c        Robert L. Parker and Philip B. Stark    Version 7/87
c
c  Solves the problem: 
c
c          min  || a.x - b ||     such that   bl <= x <= bu
c                            2
c    where  
c               x  is an unknown n-vector
c               a  is a given m by n matrix
c               b  is a given  m-vector 
c               bl is a given n-vector of lower bounds on the
c                                components of x.
c               bu is a given n-vector of upper bounds on the
c                                components of x.
c
c-----------------------------------------------------------------------
c    Input parameters:
c
c  m, n, a, b, bl, bu   see above.  It is assumed  that m < n.
c
c  If key = 0, the subroutine solves the problem from scratch.
c
c  If key > 0 the routine initializes using the user's guess about
c   which components of  x  are `active', i.e. are stricly within their
c   bounds, which are at their lower bounds, and which are at their 
c   upper bounds.  This information is supplied through the array  
c   istate.  istate(n+1) should contain the total number of components 
c   at their bounds (the `bound variables').  The absolute values of the
c   first nbound=istate(n+1) entries of  istate  are the indices
c   of these `bound' components of  x.  The sign of istate(j), j=1,...,
c   nbound, indicates whether  x(|istate(j)|) is at its upper or lower
c   bound.  istate(j) is positive if the component is at its upper
c   bound, negative if the component is at its lower bound.
c   istate(j), j=nbound+1,...,n  contain the indices of the components
c   of  x  that are active (i.e. are expected to lie strictly within 
c   their bounds).  When key > 0, the routine initially sets the active 
c   components to the averages of their upper and lower bounds: 
c   x(j)=(bl(j)+bu(j))/2, for j in the active set.  
c
c-----------------------------------------------------------------------
c    Output parameters:
c
c  x       the solution vector.
c
c  w(1)    the minimum 2-norm || a.x-b ||.
c
c  istate  vector indicating which components of  x  are active and 
c          which are at their bounds (see the previous paragraph).  
c          istate can be supplied to the routine to give it a good 
c          starting guess for the solution.
c
c  loopA   number of iterations taken in the main loop, Loop A.
c
c-----------------------------------------------------------------------
c    Working  arrays:
c
c  w      dimension n.               act      dimension m*(m+2).
c  zz     dimension m.               istate   dimension n+1.
c
c-----------------------------------------------------------------------
c  Method: active variable method along the general plan of NNLS by
c  Lawson & Hanson, "Solving Least Squares Problems," 1974.  See
c  Algorithm 23.10.  Step numbers in comment statements refer to their 
c  scheme.
c
c-----------------------------------------------------------------------
c  A number of measures are taken to enhance numerical reliability:
c
c 1. As noted by Lawson and Hanson, roundoff errors in the computation
c   of the gradient of the misfit may cause a component on the bounds
c   to appear to want to become active, yet when the component is added
c   to the active set, it moves away from the feasible region.  In this
c   case the component is not made active, the gradient of the misfit
c   with respect to a change in that component is set to zero, and the
c   program returns to the Kuhn-Tucker test.  Flag  ifrom5  is used in 
c   this test, which occurs at the end of Step 6.
c
c
c 2. When the least-squares minimizer after Step 6 is infeasible, it
c   is used in a convex interpolation with the previous solution to 
c   obtain a feasible vector.  The constant in this interpolation is
c   supposed to put at least one component of  x   on a bound. There can
c   be difficulties: 
c
c 2a. Sometimes, due to roundoff, no interpolated component ends up on 
c   a bound.  The code in Step 11 uses the flag  jj, computed in Step 8,
c   to ensure that at least the component that determined the 
c   interpolation constant  alpha  is moved to the appropriate bound.  
c   This guarantees that what Lawson and Hanson call `Loop B' is finite.
c
c 2b. The code in Step 11 also incorporates Lawson and Hanson's feature
c   that any components remaining infeasible at this stage (which must
c   be due to roundoff) are moved to their nearer bound.
c
c
c 3. If the columns of  a  passed to qr are linearly dependent, the new
c   potentially active component is not introduced: the gradient of the
c   misfit with respect to that component is set to zero, and control
c   returns to the Kuhn-Tucker test.
c
c
c 4. When some of the columns of  a  are approximately linearly 
c   dependent, we have observed cycling of active components: a 
c   component just moved to a bound desires immediately to become 
c   active again; qr allows it to become active and a different 
c   component is moved to its bound.   This component immediately wants
c   to become active, which qr allows, and the original component is
c   moved back to its bound.  We have taken two steps to avoid this 
c   problem:
c
c 4a. First, the column of the matrix  a  corresponding to the new 
c   potentially active component is passed to qr as the last column of 
c   its matrix.  This ordering tends to make a component recently moved
c   to a bound fail the test mentioned in (1), above.
c
c 4b. Second, we have incorporated a test that prohibits short cycles.
c   If the most recent successful change to the active set was to move
c   the component x(jj) to a bound, x(jj) is not permitted to reenter 
c   the solution at this stage.  This test occurs just after checking
c   the Kuhn-Tucker conditions, and uses the flag  jj, set in Step 8.
c   The flag  jj  is reset after Step 6 if Step 6 was entered from
c   Step 5 indicating that a new component has successfully entered the
c   active set. The test for resetting  jj  uses the flag  ifrom5,
c   which will not equal zero in case Step 6 was entered from Step 5.
c
c
      dimension a(m,n), b(m), x(n), bl(n), bu(n)
      dimension w(n), act(m,m+2), zz(m), istate(n+1)
      common /io/ inp,iout,iprint,ierr
c
      data eps/1.0d-11/
c%%   data eps/1.0e-05/
c
c----------------------First Executable Statement-----------------------
c
c  Step 1.  Initialize everything--active and bound sets, initial 
c   values, etc.
c
c  Initialize flags, etc.
      jj = 0
      ifrom5 = 0
      m1 = m + 1
c  Check consistency of given bounds  bl, bu.
      bdiff = 0.0
      do 1005 j=1, n
        bdiff=max(bdiff, bu(j)-bl(j))
        if (bl(j) .gt. bu(j)) then
          write(iout,*)'>>> Inconsistent bounds in BVLS. '
          stop
        endif
 1005 continue
      if (bdiff .eq. 0.0) then
       write(iout,*)'>>> No free variables in BVLS--check input bounds.'
       stop
      endif
c
c  In a fresh initialization (key = 0) bind all variables at their lower
c   bounds.  If (key != 0), use the supplied  istate  vector to
c   initialize the variables.  istate(n+1) contains the number of
c   bound variables.  The absolute values of the first 
c   nbound=istate(n+1) entries of  istate  are the indices of the bound
c   variables.  The sign of each entry determines whether the indicated
c   variable is at its upper (positive) or lower (negative) bound.
      if (key .eq. 0) then
        nbound=n
        nact=0
        do 1010 j=1, nbound
          istate(j)=-j
 1010   continue
      else
        nbound=istate(n+1)
      endif
      nact=n - nbound
      if ( nact .gt. m ) then
        write(iout,*)
     $  '>>> Too many active variables in BVLS starting solution!'
        stop
      endif
      do 1100 k=1, nbound
        j=abs(istate(k))
        if (istate(k) .lt. 0) x(j)=bl(j)
        if (istate(k) .gt. 0) x(j)=bu(j)
 1100 continue
c
c  In a warm start (key != 0) initialize the active variables to 
c   (bl+bu)/2.  This is needed in case the initial qr results in 
c   active variables out-of-bounds and Steps 8-11 get executed the 
c   first time through. 
      do 1150 k=nbound+1,n
        kk=istate(k)
        x(kk)=(bu(kk)+bl(kk))/2
 1150 continue
c
c  Compute bnorm, the norm of the data vector b, for reference.
      bsq=0.0
      do 1200 i=1, m
        bsq=bsq + b(i)**2
 1200 continue
      bnorm=sqrt(bsq)
c
c-----------------------------Main Loop---------------------------------
c
c  Initialization complete.  Begin major loop (Loop A).
      do 15000 loopA=1, 3*n
c
c  Step 2.
c  Initialize the negative gradient vector w(*).
      obj=0.0
      do 2050 j=1, n
        w(j)=0.0
 2050 continue
c
c  Compute the residual vector b-a.x , the negative gradient vector
c   w(*), and the current objective value obj = || a.x - b ||.
c   The residual vector is stored in the m+1'st column of act(*,*).
      do 2300 i=1, m
        ri=b(i)
        do 2100 j=1, n
          ri=ri - a(i,j)*x(j)
 2100   continue
        obj=obj + ri**2
        do 2200 j=1, n
          w(j)=w(j) + a(i,j)*ri
 2200   continue
        act(i,m1)=ri
 2300 continue
c
c  Converged?  Stop if the misfit << || b ||, or if all components are 
c   active (unless this is the first iteration from a warm start). 
      if (sqrt(obj) .le. bnorm*eps .or. 
     +  (loopA .gt. 1 .and. nbound .eq. 0)) then
         istate(n+1)=nbound
c  Report the norm.
         w(1)=sqrt(obj)
         return
      endif
c
c  Add the contribution of the active components back into the residual.
      do 2500 k=nbound+1, n
        j=istate(k)
        do 2400 i=1, m
          act(i,m1)=act(i,m1) + a(i,j)*x(j)
 2400   continue
 2500 continue
c
c  The first iteration in a warm start requires immediate qr.
      if (loopA .eq. 1 .and. key .ne. 0) goto 6000
c
c  Steps 3, 4.
c  Find the bound element that most wants to be active.
 3000 worst=0.0
      it=1
      do 3100 j=1, nbound
         ks=abs(istate(j))
         bad=w(ks)*sign(1, istate(j))
         if (bad .lt. worst) then
            it=j
            worst=bad
            iact=ks
         endif
 3100 continue
c
c  Test whether the Kuhn-Tucker condition is met.
      if (worst .ge. 0.0 ) then
         istate(n+1)=nbound
         w(1)=sqrt(obj)
         return
      endif
c
c  The component  x(iact)  is the one that most wants to become active.
c   If the last successful change in the active set was to move x(iact)
c   to a bound, don't let x(iact) in now: set the derivative of the 
c   misfit with respect to x(iact) to zero and return to the Kuhn-Tucker
c   test.
      if ( iact .eq. jj ) then
        w(jj)=0.0
        goto 3000
      endif
c
c  Step 5.
c  Undo the effect of the new (potentially) active variable on the 
c   residual vector.
      if (istate(it) .gt. 0) bound=bu(iact)
      if (istate(it) .lt. 0) bound=bl(iact)
      do 5103 i=1, m
        act(i,m1)=act(i,m1) + bound*a(i,iact)
 5103 continue
c
c  Set flag ifrom5, indicating that Step 6 was entered from Step 5.
c   This forms the basis of a test for instability: the gradient
c   calculation shows that x(iact) wants to join the active set; if 
c   qr puts x(iact) beyond the bound from which it came, the gradient 
c   calculation was in error and the variable should not have been 
c   introduced.
      ifrom5=istate(it)
c
c  Swap the indices (in istate) of the new active variable and the
c   rightmost bound variable; `unbind' that location by decrementing
c   nbound.
      istate(it)=istate(nbound)
      nbound=nbound - 1
      nact=nact + 1
      istate(nbound+1)=iact
c
      if (m .lt. nact) then
        write(iout,*)'>>> Too many free variables in BVLS!'
        stop
      endif
c
c  Step 6.
c  Load array  act  with the appropriate columns of  a  for qr.  For
c   added stability, reverse the column ordering so that the most
c   recent addition to the active set is in the last column.  Also 
c   copy the residual vector from act(., m1) into act(., m1+1).
 6000 do 6200 i=1, m
        act(i,m1+1)=act(i,m1)
        do 6100 k=nbound+1, n
          j=istate(k)
          act(i,nact+1-k+nbound)=a(i,j)
 6100   continue
 6200 continue
c
      call qr(m, nact, act, act(1,m1+1), zz, resq)
c
c  Test for linear dependence in qr, and for an instability that moves
c   the variable just introduced away from the feasible region 
c   (rather than into the region or all the way through it). 
c   In either case, remove the latest vector introduced from the
c   active set and adjust the residual vector accordingly.  
c   Set the gradient component (w(iact)) to zero and return to 
c   the Kuhn-Tucker test.
      if (resq .lt. 0.0 
     +   .or. (ifrom5 .gt. 0 .and. zz(nact) .gt. bu(iact))
     +   .or. (ifrom5 .lt. 0 .and. zz(nact) .lt. bl(iact))) then
         nbound=nbound + 1
         istate(nbound)=istate(nbound)*sign(1.0d0, x(iact)-bu(iact))
         nact=nact - 1
         do 6500 i=1, m
           act(i,m1)=act(i,m1) - x(iact)*a(i,iact)
 6500    continue
         ifrom5=0
         w(iact)=0.0
         goto 3000
      endif
c
c  If Step 6 was entered from Step 5 and we are here, a new variable 
c   has been successfully introduced into the active set; the last 
c   variable that was fixed at a bound is again permitted to become 
c   active.
      if ( ifrom5 .ne. 0 ) jj=0
      ifrom5=0
c
c   Step 7.  Check for strict feasibility of the new qr solution.
      do 7100 k=1, nact
        k1=k
        j=istate(k+nbound)
        if (zz(nact+1-k).lt.bl(j) .or. zz(nact+1-k).gt.bu(j)) goto 8000
 7100 continue
      do 7200 k=1, nact
        j=istate(k+nbound)
        x(j)=zz(nact+1-k)
 7200 continue
c  New iterate is feasible; back to the top.
      goto 15000
c
c  Steps 8, 9.
 8000 alpha=2.0
      alf=alpha
      do 8200 k=k1, nact
        j=istate(k+nbound)
        if (zz(nact+1-k) .gt. bu(j)) 
     +    alf=(bu(j)-x(j))/(zz(nact+1-k)-x(j))
        if (zz(nact+1-k) .lt. bl(j)) 
     +    alf=(bl(j)-x(j))/(zz(nact+1-k)-x(j))
        if (alf .lt. alpha) then
          alpha=alf
          jj=j
          sj=sign(1.0d0, zz(nact+1-k)-bl(j))
        endif
 8200 continue
c
c  Step 10
      do 10000 k=1, nact
        j=istate(k+nbound)
        x(j)=x(j) + alpha*(zz(nact+1-k)-x(j))
10000 continue
c
c  Step 11.  
c  Move the variable that determined alpha to the appropriate bound.  
c   (jj is its index; sj is + if zz(jj)> bu(jj), - if zz(jj)<bl(jj) ).
c   If any other component of  x  is infeasible at this stage, it must
c   be due to roundoff.  Bind every infeasible component and every
c   component at a bound to the appropriate bound.  Correct the
c   residual vector for any variables moved to bounds.  Since at least
c   one variable is removed from the active set in this step, Loop B 
c   (Steps 6-11) terminates after at most  nact  steps.
      noldb=nbound
      do 11200 k=1, nact
        j=istate(k+noldb)
        if (((bu(j)-x(j)) .le. 0.0) .or. 
     +    (j .eq. jj .and. sj .gt. 0.0)) then
c  Move x(j) to its upper bound.
          x(j)=bu(j)
          istate(k+noldb)=istate(nbound+1)
          istate(nbound+1)=j
          nbound=nbound+1
          do 11100 i=1, m
             act(i,m1)=act(i,m1) - bu(j)*a(i,j)
11100     continue
        else if (((x(j)-bl(j)) .le. 0.0) .or.
     +    (j .eq. jj .and. sj .lt. 0.0)) then
c  Move x(j) to its lower bound.
          x(j)=bl(j)
          istate(k+noldb)=istate(nbound+1)
          istate(nbound+1)=-j
          nbound=nbound+1
          do 11150 i=1, m
             act(i,m1)=act(i,m1) - bl(j)*a(i,j)
11150     continue
        endif
11200 continue
      nact=n - nbound
c
c  If there are still active variables left repeat the qr; if not,
c    go back to the top.
      if (nact .gt. 0 ) goto 6000
c
15000 continue
c
      write(iout,*)'>>> BVLS fails to converge! '
      stop
      end
c_____________________________________________________________________

      subroutine qr(m, n, a, b, x, resq)
      implicit double precision (a-h, o-z)
c$$$$ calls no other routines
c  Relies on FORTRAN77 do-loop conventions!
c  Solves over-determined least-squares problem  ax ~ b
c  where  a  is an  m by n  matrix,  b  is an m-vector .
c  resq  is the sum of squared residuals of optimal solution.  Also used
c  to signal error conditions - if -2 , system is underdetermined,  if
c  -1,  system is singular.
c  Method - successive Householder rotations.  See Lawson & Hanson - 
c  Solving Least Squares Problems (1974).
c  Routine will also work when m=n.
c*****   CAUTION -  a and b  are overwritten by this routine.
      dimension a(m,n),b(m),x(n)
c  Even if the real variables are used elsewhere, retain the statement 
c  below for improved precision in the solution.
      double precision sum, dot
      common /io/ inp,iout,iprint,ierr
c
      resq=-2.0
      if (m .lt. n) return
      resq=-1.0
c   Loop ending on 1800 rotates  a  into upper triangular form.
      do 1800 j=1, n
c   Find constants for rotation and diagonal entry.
        sq=0.0
        do 1100 i=j, m
          sq=a(i,j)**2 + sq
 1100   continue
        if (sq .eq. 0.0) return
        qv1=-sign(sqrt(sq), a(j,j))
        u1=a(j,j) - qv1
        a(j,j)=qv1
        j1=j + 1
c  Rotate remaining columns of sub-matrix.
        do 1400 jj=j1, n
          dot=u1*a(j,jj)
          do 1200 i=j1, m
            dot=a(i,jj)*a(i,j) + dot
 1200     continue
          const=dot/abs(qv1*u1)
          do 1300 i=j1, m
            a(i,jj)=a(i,jj) - const*a(i,j)
 1300     continue
          a(j,jj)=a(j,jj) - const*u1
 1400   continue
c  Rotate  b  vector.
        dot=u1*b(j)
        do 1600 i=j1, m
          dot=b(i)*a(i,j) + dot
 1600   continue
        const=dot/abs(qv1*u1)
        b(j)=b(j) - const*u1
        do 1700 i=j1, m
          b(i)=b(i) - const*a(i,j)
 1700   continue
 1800 continue
c  Solve triangular system by back-substitution.
      do 2200 ii=1, n
        i=n-ii+1
        sum=b(i)
        do 2100 j=i+1, n
          sum=sum - a(i,j)*x(j)
 2100   continue
        if (a(i,i).eq. 0.0) return
         x(i)=sum/a(i,i)
 2200 continue
      resq=m - n
* Since  only the sign of  resq is used in bvls, following lines
* are sterilized
cc  Find residual in overdetermined case.
c       resq=0.0
c       do 2300 i=n+1, m
c         resq=b(i)**2 + resq
c  2300 continue
      return
      end

c_______________________________________________________________________
 
      function fit(nrefin)
      implicit double precision (a-h, o-z)

c     nrefin is number of refinement iterations (1 usually adaquate)
 
c  does rho+ fit with for data included according to drh, dph >0
c
c     iband is the band to offset (2 to 4)
c     Note: 1 is the fixed band and is NEVER offset
c
      parameter (inmx=200)
      parameter (mxf=101,mxxi=450, mxrw=2*mxf)
      common // a(mxrw*mxxi),d(mxrw),x(mxxi),ma,na
      common /bounds/ bl(mxxi),bu(mxxi)
      common /work/ w(mxxi),zz(mxrw),act(mxrw,mxrw+2),istt(mxxi+1)
      common /shared/ pi,twopi,amu0,radian
      common /observ/ freq(mxf),rho(mxf),drh(mxf),phi(mxf),dph(mxf)
     $  ,erh(mxf),eph(mxf),dwtrh(mxf),dwtph(mxf),sgr(mxf),nfrq
      common /final/ respon(2,mxf),bfresp(2,mxf),bound(2,2*mxf)
      common /discrt/ xi(mxxi),nxi
      common /optiml/ wtmax,fitmin,prob
      common /prolix/ loquor
      common /choice/ a0, state
      common /io/ inp,iout,iprint,ierr

      character*10 state

      common /bopram/ iband,nbands,nstart(4),nend(4),gain(4)
      logical logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs(4)
      character*3 cecho
      common /flags/ logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs,cecho


      call sample()
      call builda(ma, nfrq, a, d)
      call bvls(0, ma,nxi,a,d, bl,bu, x, w,act,zz,istt, loopA)
      fit=w(1)**2
      if(loquor.ge.2) write(iout,*)'Initial fit= ',fit
      if (loquor .gt. 3) call report()
      do 100 i=1,nrefin
        fit1=fit
        call refine()
        call builda(ma, nfrq, a, d)
        call bvls(0, ma,nxi,a,d, bl,bu, x, w,act,zz,istt, loopA)
        fit=w(1)**2
        if(loquor.ge.2)write(iout,*)'Refinement phase ',i,': fit= ',fit
        if (fit .gt. fit1) then
           write(iout,'(2a,i2,a,i2,a)')'  ',
     $'>>> Warning: Refinement phase',i,' of ',nrefin,' failed to improv
     $e fit'
           write(iout,'(a,g17.9,a,g17.9)')
     $     'Chisq before: ',fit1,'  Chisq after: ',fit
        endif
 100  continue
      
      if(flg3) write(iout, '(2a)')'Surface condition: ',state

      return
      end

c_______________________________________________________________________

      function fitgn(tgain)
      implicit double precision (a-h, o-z)

c  does rho+ fit with data in a band set in common /bopram/
c
c     iband is the band to offset (2 to 4)
c     Note: 1 is the fixed band and is NEVER offset
c  
      parameter (inmx=200)
      parameter (mxf=101,mxxi=450, mxrw=2*mxf)
      common // a(mxrw*mxxi),d(mxrw),x(mxxi),ma,na
      common /bounds/ bl(mxxi),bu(mxxi)
      common /work/ w(mxxi),zz(mxrw),act(mxrw,mxrw+2),istt(mxxi+1)
      common /shared/ pi,twopi,amu0,radian
      common /observ/ freq(mxf),rho(mxf),drh(mxf),phi(mxf),dph(mxf)
     $  ,erh(mxf),eph(mxf),dwtrh(mxf),dwtph(mxf),sgr(mxf),nfrq
      common /final/ respon(2,mxf),bfresp(2,mxf),bound(2,2*mxf)
      common /discrt/ xi(mxxi),nxi
      common /optiml/ wtmax,fitmin,prob
      common /prolix/ loquor
      common /io/ inp,iout,iprint,ierr
      common /choice/ a0, state

      character*10 state

      common /bopram/ iband,nbands,nstart(4),nend(4),gain(4)
      logical logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs(4)
      character*3 cecho
      common /flags/ logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs,cecho


      common /save/savrho(mxf),savdrh(mxf),savphi(mxf),savdph(mxf)

      dimension tmprho(mxf),tmpdrh(mxf)

      if(iband.eq.1)stop 'illegal iband=1 in function fitgn'

c save rhos in chosen band and then offset them
c also rescale errors so that % log error stays the same

      do 10 i=nstart(iband),nend(iband)
         tmprho(i)=rho(i)
         tmpdrh(i)=drh(i)
         g=10.**tgain
         rho(i)=rho(i)*g
         drh(i)=drh(i)*g
10    continue

      fitmin=fit(1)
c     fitmin=fit(0)

c if using only phase fit to judge offset, calculate phase chisq for
c   the included phase data only
      if(phjudg)then
        call forwrd(nfrq, freq, respon)
        fitph=0.
        do 50 i=1,nfrq
          if(savdph(i).gt.0.)then
             fitph=fitph+((respon(2,i)-savphi(i))/savdph(i))**2
          endif
50      continue
        fitgn=fitph
      else
c otherwise use chisq of all included data
        fitgn=fitmin
      endif

c restore rhos and errors in offset band
      do 100 i=nstart(iband),nend(iband)
         rho(i)=tmprho(i)
         drh(i)=tmpdrh(i)
100   continue

      return
      end

c_______________________________________________________________________

      subroutine fndofs()

      implicit double precision (a-h, o-z)
c
      parameter (inmx=200)
      character *80 input(inmx)
      common /store/ input
      common /ndict/ iecho,nin,memory,istate(inmx)
c
      parameter (mxf=101,mxxi=450, mxrw=2*mxf)
      common // a(mxrw*mxxi),d(mxrw),x(mxxi),ma,na
      common /bounds/ bl(mxxi),bu(mxxi)
      common /work/ w(mxxi),zz(mxrw),act(mxrw,mxrw+2),istt(mxxi+1)
      common /shared/ pi,twopi,amu0,radian
      common /observ/ freq(mxf),rho(mxf),drh(mxf),phi(mxf),dph(mxf)
     $  ,erh(mxf),eph(mxf),dwtrh(mxf),dwtph(mxf),sgr(mxf),nfrq
      common /final/ respon(2,mxf),bfresp(2,mxf),bound(2,2*mxf)
      common /discrt/ xi(mxxi),nxi
      common /optiml/ wtmax,fitmin,prob
      common /prolix/ loquor
      common /io/ inp,iout,iprint,ierr
      common /save/savrho(mxf),savdrh(mxf),savphi(mxf),savdph(mxf)

      character*10 state
      common /choice/ a0,state
      common /bopram/ iband,nbands,nstart(4),nend(4),gain(4)

      logical logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs(4)
      character*3 cecho
      common /flags/ logdat,phjudg,flg1,flg2,flg3,rhoflg,offst,
     $               usephi,phfit,fixgn,incofs,cecho

      external fitgn,mnbrak,brent

c  Exclude phase if not used in gain calculation

         if(.not.usephi)then
            do 15 i=1,nfrq
               dph(i)=0.
 15         continue
         endif

c  Save original data and error scales

C        do 10 i=1,nfrq
C           savrho(i)=rho(i)
C           savdrh(i)=drh(i)
C           savphi(i)=phi(i)
C           savdph(i)=dph(i)
C10      continue
 

         if(usephi)then
           write(ierr,*)' '
           write(ierr,*)'Rho and phase used in band gain calculation'
           write(iout,*)
           write(iout,*)'Rho and phase used in band gain calculation'
         else
           write(ierr,*)' '
           write(ierr,*)'Only rho used in band gain calculation'
           write(iout,*)'Only rho used in band gain calculation'
         endif

         if(phjudg)then
           write(ierr,*)'Gains judged using fit to phase only'
           write(iout,*)'Gains judged using fit to phase only'
         elseif(usephi)then
           write(ierr,*)'Gains judged using fit to rho & phase'
           write(iout,*)'Gains judged using fit to rho & phase'
         endif

         if(nbands.gt.2)then
           if(rhoflg )then
            write(ierr,*)'Using adjusted rho from earlier fit bands'
            write(iout,*)'Using adjusted rho from earlier fit bands'
           else
            write(ierr,*)'Using rho in fixed, plus band being worked on'
            write(iout,*)'Using rho in fixed, plus band being worked on'
           endif
         endif

         if( (.not.usephi) .and. phfit)then
           write(ierr,*)'Phase used in fits after gain adjustment'
           write(iout,*)'Phase used in fits after gain adjustment'
         endif

         tol=1e-3
         if(loquor .ge. 3)
     &     write(iout,'(a,f10.9)')'tolerance for brent= ',tol

c  refine estimates of band offsets

         do 60 k=2,nbands

           if(incofs(k))then

c  modify errors in rho to include fixed band(s) and 2 through kth band

             if(rhoflg)then
               do 20 i=1,nfrq
                  drh(i)=0.
 20            continue
               do 22 i=nstart(1),nend(1)
                  drh(i)=savdrh(i)
 22            continue
               if(k.ge.2 .or. .not. incofs(k))then
                  do 24 i=nstart(2),nend(2)
                     drh(i)=savdrh(i)
 24               continue
               endif
               if(k.ge.3 .or. .not.incofs(k))then
                  do 26 i=nstart(3),nend(3)
                     drh(i)=savdrh(i)
 26               continue
               endif
               if(k.eq.4 .or. .not.incofs(k))then
                  do 28 i=nstart(4),nend(4)
                     drh(i)=savdrh(i)
 28               continue
               endif

             else

c  modify errors in rho to include only fixed and kth bands

               do 30 i=1,nfrq
                  drh(i)=0.
 30            continue
               do 32 i=nstart(1),nend(1)
                  drh(i)=savdrh(i)
 32            continue
               do 34 i=nstart(k),nend(k)
                  drh(i)=savdrh(i)
 34            continue
               if(.not.incofs(k))then
                  do 36 i=nstart(k),nend(k)
                     drh(i)=savdrh(i)
 36               continue
               endif
             endif

             if(loquor .ge. 3)then
               write(iout,*)
     &         'savrho(j),rho(j),savdrh(j),drh(j),phi(j),dph(j):'
               do 51 i=1,nfrq
                 write(iout,'(i3,6g12.4)')
     &           i,savrho(i),rho(i),savdrh(i),drh(i),phi(i),dph(i)
 51            continue
             endif

c  suppress printing of sampling notices by sample and
c  refine after first iteration

             flg1=.true.
             iband=k

c  find starting triple bracket for gain 

             x1=  0.
             x2= -0.1
             call mnbrak(x1,x2,x3,f1,f2,f3,fitgn)

c  use Brent algorithm to find band gain that minimizes chisq

             fitmp=brent(x1,x2,x3,fitgn,tol,gain(k))
c            write(iout,'(a,g12.4)')'fitmp= ',fitmp
c            write(ierr,'(4g12.4)')(10.**gain(i),i=1,4)

c  apply new gain permanently to kth band
             do 50 i=nstart(k),nend(k)
                savrho(i)=savrho(i)*(10.**gain(k))
                savdrh(i)=savdrh(i)*(10.**gain(k))
                rho(i)=savrho(i)
                drh(i)=savdrh(i)
 50          continue
           endif
 60      continue

         write(iout,*)
         if(logdat)then
          write(iout,*)
     &      'Calculated optimum log10 gains added to log10(rho) data:'
          write(iout,'(a6,i3,a2,i3,1x,f7.4,a6)')' Band ',nstart(1),' :',
     &                              nend(1),gain(1),' Fixed'
          do 65 k=2,nbands
            write(iout,'(a6,i3,a2,i3,1x,f7.4)')' Band ',nstart(k),' :',
     &                              nend(k),gain(k)
 65       continue
         else
          write(iout,*)
     &      'Calculated optimum gains multiplied times rho data:'
          write(iout,'(a6,i3,a2,i3,1x,f7.4,a6)')' Band ',nstart(1),' :',
     &                              nend(1),10.**gain(1),' Fixed'
          do 66 k=2,nbands
            write(iout,'(a6,i3,a2,i3,1x,f7.4)')' Band ',nstart(k),' :',
     &                              nend(k),10.**gain(k)
 66       continue
         endif

c restore rho to its saved state
         do 68 i=1,nfrq
            rho(i)=savrho(i)
            drh(i)=savdrh(i)
 68      continue

c  Restore phase if it was not used in gain calculation

         if(.not.usephi)then
            do 75 i=1,nfrq
               dph(i)=savdph(i)
 75         continue
         endif

         return
         end
c_______________________________________________________________________
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c The following two routines "brent" and "mnbrak"
c   are from:
c
c   Numerical Recipies: The Art of Scientific Computing
c   Press, Flannery, Teulosky and Vetterling
c   Cambridge University Press
c
c These routines are copyrighted and their use presumes that you
c   have purchased the book and/or the diskettes that go with it.
c   (A VERY worthwile book.)
c
c These routines are used only in the section that computes the
c   best offset for a tear in the apparent resistivity. If
c   you do not have a license to use these routines, you should
c   set numrec=.false. in the main program. This will disable the
c   calculation of optimum tear offset part of the code. The rest of
c   the code, including manual setting of gains in sub-bands is
c   still enabled.
c
c If you do not have permission to use these routines, you can use
c   replacements. They simply find the minimum of a function of a
c   single variable in an efficient way. If you do so and they work well,
c   please pass the routines on to us (JRB and RLP.
c_______________________________________________________________________

       function brent(ax,bx,cx,f,tol,xmin)

       IMPLICIT DOUBLE PRECISION (a-h, o-z)

c inputs:
c     real ax, bx, cx, f(x), tol
c outputs:
c     real brent, xmin
c
c     Brent's algorithm for finding the minimum of a function
c        from Numerical Recipies
c
c      ax < bx < cx abcissas of function f(x) such that ax and cx
c        bracket the minimum of f.
c
c      The minimum is brent = f(xmin), where f(x) is a user-supplied
c         external function
c
c      WARNING!!! Make sure that the function f (probably will
c        have a different name) is declared EXTERNAL by the
c        calling program. Otherwise this function will bomb.

c      real cgold, zeps
       parameter(cgold=.3819660, zeps=1.0d-10)

       integer itmax
       parameter(itmax=100)

c      real a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
       common /io/ inp,iout,iprint,ierr
       common /prolix/ loquor

       integer iter

       external f
c
       a=min(ax,cx)
       b=max(ax,cx)
       v=bx
       w=v
       x=v
       e=0.
       fx=f(x)
       fv=fx
       fw=fx

c  Main loop

       do 100 iter=1,itmax
          xm=0.5*(a+b)
          tol1=tol*abs(x)+zeps
          tol2=2.*tol1

c  Check to see if you have converged

          if(abs(x-xm).le.(tol2 - 0.5*(b-a)))then
             xmin=x
             brent=fx
             if(loquor.ge.1) then
C              write(iout,*)'Brent converged at ',iter,' iterations'
               write(iout,*)
     &            'Gain search converged at ',iter,' iterations'
             endif
             return
          endif

          if(abs(e).gt.tol1)then
             r=(x-w)*(fx-fv)
             q=(x-v)*(fx-fw)
             p=(x-v)*q-(x-w)*r
             q=2.*(q-r)
             if(q.gt.0.) p=-p
             q=abs(q)
             etemp=e
             e=d

c  If a paraolic step NOT OK go to a golden section step

             if(abs(p).ge.abs(0.5*q*etemp) .or. p.le.q*(a-x) .or.
     &          p.ge.q*(b-x)) go to 1

c  Take parabolic step

             d=p/q
             u=x+d
             if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
             go to 2
          endif
  
c  or take golden section step

 1        if(x.ge.xm)then
             e=a-x
          else
             e=b-x
          endif
          d=cgold*e

 2        if(abs(d).ge.tol1)then
             u=x+d
          else
             u=x+sign(tol1,d)
          endif
          fu=f(u)
          if(fu.le.fx)then
             if(u.ge.x)then
                a=x
             else
                b=x
             endif
             v=w
             fv=fw
             w=x
             fw=fx
             x=u
             fx=fu
          else
             if(u.lt.x)then
                a=u
             else
                b=u
             endif
             if(fu.le.fw .or. w.eq.x)then
                v=w
                fv=fw
                w=u
                fw=fu
             elseif(fu.le.fv .or. v.eq.x .or. v.eq.w)then
                v=u
                fv=fu
             endif
          endif
          
 100   continue

c  End of main loop

c  If you get here:

       write(iout,*)'Warning Brent minimum-finder exceeded ',
     &           itmax,' iterations'

       xmin=x
       brent=fx

       return
       end

c_______________________________________________________________________

       subroutine mnbrak(ax,bx,cx,fa,fb,fc,func)
       implicit double precision (a-h, o-z)

c  Brackets a minimum of func(x). Used to start function brent, if
c    a good start is unknown. 

c  From Numerical Recipies by Press et al.

c Input:
c  ax,bx starting guesses. (Should be as close as feasible, but need
c    not bracket minimum.)

c Output:
c  ax,bx,cx, where ax and cx bracket the minimum and fa=func(ax),
c    fb=func(bx) etc.
c
c  WARNING!!! Make sure that the function func (probably will
c    have a different name) is declared EXTERNAL by the
c    calling program. Otherwise this subroutine will bomb.
c
       external func
       parameter(gold=1.618034, glimit=100, tiny=1.d-20)

       fa=func(ax)
       fb=func(bx)

c switch roles of ax and bx if needed to ensure downhill search

       if(fb.gt.fa)then
          dum=ax
          ax=bx
          bx=dum
          dum=fb
          fb=fa
          fa=dum
       endif

c first guess, linear 'golden section' extrapolation

       cx=bx+gold*(bx-ax)
       fc=func(cx)

c start of main loop

c if we have not yet passed minimum, do parabolic interpolation/
c   extrapolation to estimated minimum point ux.

 1     if(fb.ge.fc)then
          r=(bx-ax)*(fb-fc)
          q=(bx-cx)*(fb-fa)
          ux=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),tiny),q-r))

c (the presence of tiny presvents accidental division by zero,
c   should the initial guesses accidently be precisely equal
c   distances on either side of the parabola's minimum.)

c This will be largest distance that one can go if ux > cx.

          ulim=bx+glimit*(cx-bx)

c check to see whether ux is between bx and cx


          if((bx-ux)*(ux-cx).gt.0.)then
             fu=func(ux)

c test to see if ux is smaller than bx or cx, if so, then bx,ux,cx
c   is the desired bracket

             if(fu.lt.fc)then
                ax=bx
                fa=fb
                bx=ux
                fb=fu
                return

c if that doesn't work see if there is a minimum between ux and ax.
c   If so, then ax,bx,ux is the desired bracket.

             elseif(fu.gt.fb)then
                cx=ux
                fc=fu
                return
             endif

c if you get here, parabolic interpolation did not lead to a
c   ux that did the job. Try a larger "golden" step.

             ux=cx+gold*(cx-bx)
             fu=func(ux)

c if ux > cx see if there is a minimum between cx and ulim.

          elseif((cx-ux)*(ux-ulim).ge.0.)then

             fu=func(ux)
             if(fu.lt.fc)then
c you are still going downhill try another "golden step"
                bx=cx
                cx=ux
                ux=cx+gold*(cx-bx)
                fb=fc
                fc=fu
                fu=func(ux)
             endif

c make sure parabolic extrapolation does not get too far away
c   from present guesses.

          elseif((ux-ulim)*(ulim-cx).ge.0.)then
             ux=ulim
             fu=func(ux)

c back to successively amplified "golden" steps

          else
             ux=cx+gold*(cx-bx)
             fu=func(ux)
          endif

c eliminate oldest point and continue

          ax=bx
          bx=cx
          cx=ux
          fa=fb
          fb=fc
          fc=fu
          go to 1
c end of main loop

       endif

       return
       end
c_______________________________________________________________________
c
c End of routines from Numerical Recipies
c
c_______________________________________________________________________

