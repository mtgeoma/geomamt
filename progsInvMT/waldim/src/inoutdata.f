C   ******************************************************************
C     Subroutine ReadParams 
C     Reads program parameters 
C     a) from file param.cfg, or 
C     b) asks the user and creates param.cfg file with these parameters

      SUBROUTINE ReadParams(file11,units,th,thq,eif,uef,ep,
     +                     avg,pm,pma,bpd,out,statf,errcfg)
      
      implicit double precision (a-h,o-z)
      Character file11*30,Info*100,eif*1,uef*1
      Character Units*1,avg*1,statf*1,out*25

      open (unit=1,file='param.cfg',status='old',err=11)

C     Reading parameters from old file 'param.cfg':
      read(1,'(a)',err= 21,end=21) info
      read(1,*,err= 21,end=21) file11

      open(unit=11,file=file11,status='old',err=15)
      close (11)
      goto 16
    
   15    err11=1
         errcfg=1

   16 read(1,'(a)',err= 21,end=21) info
      read(1,*,err= 21,end=21) Units

      if ((Units.ne.'').and.(Units.ne.'M').and.(Units.ne.'m').and.
     +(Units.ne.'F').and.(Units.ne.'f').and.(Units.ne.'Z').and.
     +(Units.ne.'z')) then
         errunits=1
         errcfg=1
      end if

      read(1,'(a)',err= 21,end=21) info
      read(1,*,err= 21,end=21) th

      if ((th.lt.0).or.(th.gt.1)) then
         errth=1
         errcfg=1
      end if

      read(1,'(a)',err= 21,end=21) info
      read(1,*,err= 21,end=21) thq

      if ((thq.lt.0).or.(thq.gt.1)) then
         errthq=1
         errcfg=1
      end if
      read(1,'(a)',err= 21,end=21) info
      read(1,*,err= 21,end=21) eif
      
      if (((eif.ne.'y').and.(eif.ne.'Y').
     +and.(eif.ne.'n').and.(eif.ne.'N'))) then 
         erreif=1
         errcfg=1
      end if

      if ((eif.eq.'y').or.(eif.eq.'Y')) then
         read(1,'(a)',err= 21,end=21) info
         read(1,*,err= 21) uef
      end if

        if ((eif.eq.'n').or.(eif.eq.'N')) then
           read(1,'(a)',err= 21,end=21) info
           read(1,'(a)',err= 21,end=21) info
           uef='n'
        end if

        if (((uef.ne.'y').and.(uef.ne.'Y').
     +     and.(uef.ne.'n').and.(uef.ne.'N'))) then 
           erruef=1
           errcfg=1
        end if
      
        if ((uef.eq.'n').or.(uef.eq.'N')) then 
           read(1,'(a)',err= 21,end=21) info
           read(1,*,err= 17,end=21) ep
           goto 18   
   17      errep=1
           errcfg=1
        else
           read(1,'(a)',err= 21,end=21) info
           read(1,'(a)',err= 21,end=21) info
           ep=0
        end if
      
   18 read(1,'(a)',err= 21,end=21) info
      read(1,*,err= 21,end=21) avg
                                  
      if ((avg.eq.'y').or.(avg.eq.'Y')) then
         read(1,'(a)',err= 21,end=21) info 
         read(1,*,err= 21,end=21) PM
         read(1,'(a)',err= 21) info
         read(1,*,err= 21,end=21) PMA
         read(1,'(a)',err= 21,end=21) info
         read(1,*,err= 21,end=21) bpd
      end if 


        if (pma.lt.pm) then
           errper=1
           errcfg=1
        end if

        if ((avg.eq.'n').or.(avg.eq.'N')) then
           read(1,'(a)',err= 21,end=21) info
           read(1,'(a)',err= 21,end=21) info
           read(1,'(a)',err= 21,end=21) info
           read(1,'(a)',err= 21,end=21) info
           read(1,'(a)',err= 21,end=21) info
           read(1,'(a)',err= 21,end=21) info  
           
           PM=0
           PMA=0
           bpd=0
        end if  

        if (((avg.ne.'y').and.(avg.ne.'Y').
     +     and.(avg.ne.'n').and.(avg.ne.'N'))) then 
           erravg=1
           errcfg=1
      end if

      read(1,'(a)',err= 21,end=21) info     
      read(1,'(a)',err= 21,end=21) out
      read(1,'(a)',err= 21,end=21) info
      read(1,'(a)',err= 21,end=21) statf

      if(((statf.ne.'y').and.(statf.ne.'Y').
     +    and.(statf.ne.'n').and.(statf.ne.'N'))) then 
          errsta=1
          errcfg=1
      end if
     
      
C     Reporting errors when reading param.cfg:
   
      goto 22

   21 errcfg=1

   22 if(err11.eq.1) then
        write(6,*)
        write(6,*) 'Error reading sites list, file ',file11,
     +' does not exist !!!!!!'
      endif

      if(errunits.eq.1) then
         write(6,*) 'Error reading units (it should be M, F or Z)!!!'
      end if

      if(errth.eq.1) then
        write(6,*)
        write(6,*) 'Error in threshold value (it should be in the range 
     +[0 - 1] !!!'
      end if

      if(errthq.eq.1) then
        write(6,*)
        write(6,*) 'Error in Q threshold value (it should be in the range 
     +[0 - 1] !!!'
      end if

      if (erreif.eq.1) then
         write(6,*)
         write(6,*) 'Not understood if EDI file contains errors 
     +(it should be Y or N) !!!'
      end if

      if(erruef.eq.1) then
        write(6,*)
        write(6,*) 'Not understood if errors in EDI files should be used 
     +(it should be Y or N) !!!' 
      end if
 
      if(errep.eq.1) then
        write(6,*)
        write(6,*) 'Error percentage should be a number !!!'
      end if
 
  
      if(erravg.eq.1) then
        write(6,*)
        write(6,*) 'Not understood if results should be averaged in bands
     +(it should be Y or N) !!!' 
      end if  
 
      if(errsta.eq.1) then
        write(6,*)
        write(6,*) 'Not understood if file with all param. and errors '
        write(6,*) 'should be created (it should be Y or N) !!!' 
      end if   
  
      if(errper.eq.1) then
        write(6,*)
        write(6,*) 'Maximum period lower than minimum period !!! '
      end if   

      if (errcfg.eq.1) then
        write(6,*) 
        write(6,*) 'Error reading param.cfg'
        write(6,*) 'Program will close'
        write(6,*) 'Revise or delete file and start again'
        write(6,*)
      end if


      goto 12


C  If file param.cfg does not exist:
C  Asking for parameters and creating file param.cfg:

   11 write(6,*)'Configuration file "param.cfg" not found.'
      write(6,*)'Creating file "param.cfg" to write parameters into.'
      write(6,*)
      open (unit=1,file='param.cfg',status='new')
               
      write(6,*) 'Name of file with edifiles names and coordinates'
      write(6,*) '[Default: list.dat]'
   24 read(5,'(a)') file11
      if (file11.eq.'') then
         file11='list.dat'
      end if
      
        open(unit=11,file=file11,status='old',err=25)
        close (11)
        goto 26

   25   write(6,*) 'File ',file11, ' does not exist'
        write(6,*) 'Please, type new file name [Default: list.dat]'
        goto 24
        
   26   write(1,'(a)') '  File with edifiles names and coordinates 
     +[Default: list.dat]='
        write(1,'(a)') file11
     
        write(6,*) 'Units for impedances in edi files?'
        write(6,'(a)') 'M=m/s, F=km/s (=(mV/km/nT) field units), Z=Ohm  
     +[Default: F (km/s)]'
   23 read(5,'(a)') Units 
            
        if (Units.eq.'') then
           Units='F'
        end if
        write(1,99) Units

        if ((Units.ne.'M').and.(Units.ne.'m').and.
     +     (Units.ne.'F').and.(Units.ne.'f').and.(Units.ne.'Z').and.
     +     (Units.ne.'z')) then
           write(6,*) 'Not understood type of units'
           write(6,*) 'Please, type units again'
           goto 23
        end if        
        
        
  271 write(6,*) 'Threshold value for invariants I3 to I7? 
     +[range 0 - 1] [Default: 0.15]'
   27   read(5,'(a)') info
           if (info.eq.'') then
             info='0.15'
           end if
                   read(info,*,err=271) th
           
        if ((th.lt.0).or.(th.gt.1)) then
         write (6,*) 'Threshold not within the range [0 - 1]'
         write(6,*) 'Please, type threshold value again [Default: 0.15]'
         goto 27
        end if


  281  write(6,*) 'Q Threshold value? [range 0 - 1] [Default: 0.1]'
   28   read(5,'(a)') info
        if (info.eq.'') then
           info='0.1'
        end if
        read(info,*,err=281) thq

        if ((thq.lt.0).or.(thq.gt.1)) then
        write (6,*) 'Q threshold not within the range [0 - 1]'
        write(6,*) 'Please, type Q threshold value again [Default: 0.1]'
        goto 28
        end if
        write(6,*) 'EDI files contain errors (Y/N)? [Default: Y]'
   29   read(5,'(a)') eif
        if (eif.eq.'') then
                   eif='Y'
                end if
      
        if (((eif.ne.'y').and.(eif.ne.'Y').
     +     and.(eif.ne.'n').and.(eif.ne.'N'))) then 
           write(6,*) 'Not understood, please type Y or N [Default: Y]'
             goto 29
        end if

        write(1,100) th,thq,eif
      
        if((eif.eq.'y').or.(eif.eq.'Y')) then
             write(6,*) 'Using errors in EDI files (Y/N)? [Default: Y]'
   30          read(5,'(a)') uef
             if (uef.eq.'') then
               uef='Y'
             end if

          if(((uef.ne.'y').and.(uef.ne.'Y').
     +    and.(uef.ne.'n').and.(uef.ne.'N'))) then 

            write(6,*) 'Not understood, please type Y or N [Default: Y]'
            goto 30
          end if
             write(1,101) uef
         else
             uef='n'
                   write(1,'(a)') '  Using errors in files (Y/N):'
             write(1,'(a)') 'NA'
         end if
        
   31          if((uef.eq.'n').or.(uef.eq.'N')) then 
                 write(6,*) 'Error percentage? [Default: 5 (5%)]'
                 read(5,'(a)') info
                 if (info.eq.'') then
                 info='5.0'
                 end if
                 read(info,*,err=31) ep
             write(1,102) ep
             else
             write(1,'(a)') '  Error percentage='
             write(1,'(a)') 'NA'
             ep=0
             end if
        write(6,*) 'Do you want to average in bands (Y/N)? [Default: Y]'
   32        read(5,'(a)') avg
             if (avg.eq.'') then
             
             avg='Y'
             end if

      if(((avg.ne.'y').and.(avg.ne.'Y').
     +    and.(avg.ne.'n').and.(avg.ne.'N'))) then 

          write(6,*) 'Not understood, please type Y or N [Default: Y]'
            goto 32
          end if

          write(1,103) avg
          
          if((avg.eq.'y').or.(avg.eq.'Y')) then
  35       write(6,*) 'Minimum period to average? [Default: 0.001]'
  33            read(5,'(a)') info
          if (info.eq.'') then
             info='0.001'
            end if
             read(info,*,err=33) pm

            write(6,*) 'Maximum period to average? [Default: 10000]'
  34          read(5,'(a)') info
             if (info.eq.'') then
             info='10000'
             end if
             read(info,*,err=34) pma
             

      if (pma.lt.pm) then

        write(6,*) 'Maximum period lower than minimum !!!'
        write(6,*) 'Please type minimum and maximum periods again'
      goto 35
      end if
             write(1,104) pm,pma
      

          write(6,*) 'Number of bands per decade? [Default: 1]'
   36     read(5,'(a)') info
          if (info.eq.'') then
             info='1'
          end if
          read(info,*,err=36) bpd
               write(1,105) bpd
               else
               write(1,'(a)') '  Minimum period to average='
               write(1,'(a)') 'NA'
               write(1,'(a)') '  Maximum period to average='
               write(1,'(a)') 'NA'
               write(1,'(a)') '  Number of bands per decade='
               write(1,'(a)') 'NA'
          end if
               
               
        write(6,*) 'Type root name for output files: [Default: OUT]'
   37   read(5,'(a)') out
        if (out.eq.'') then
        out='OUT'
        end if
        write(1,106) out

  371 write(6,*) 'Writing file with all parameters and errors (Y/N)?
     + [Default: Y]'
   38   read(5,'(a)') statf
       if (statf.eq.'') then
       statf='Y'
       end if
       
       if(((statf.ne.'y').and.(statf.ne.'Y').
     +    and.(statf.ne.'n').and.(statf.ne.'N'))) then 
          write(6,*) 'Not understood, please type Y or N [Default: Y]'
            goto 38
      end if
      write(1,107) statf


   99 format('  Units for impedances in edi files (M=m/s, F=km/s 
     +(=(mV/km/nT) field units), Z=Ohm) [Default: F(km/s)]=',/,a1)
  100 format('  Threshold value for invariants I3 to I7 
     +[range 0 - 1] [Default: 0.15]=',/,f5.3,/,
     +'  Q threshold value [range 0 - 1] [Default: 0.10]=',/,f5.3,/,
     +'  EDI files contain errors (Y/N):',/,a1)
  101 format('  Using errors in EDI files (Y/N):',/,a1)
  102 format('  Error percentage=',/,e9.3)
  103 format('  Average in bands (Y/N):',/,a1)
  104 format('  Minimum period to average=',/,e9.3,/,
     +'  Maximum period to average=',/,e9.3)
  105 format('  Number of bands per decade=',/,f4.1)
  106 format('  Root name for output files [Default: OUT]=',/,a25)
  107 format('  Writing file with all params. and errors (Y/N):',/,a1)
  
  
   12 close(1)
      return  !end of ReadParams
      
      end

C   ****************************************************************************
C  Subroutine readtbl. It reads info about each site, in any format.

      subroutine readtbl(str,p1,p2,sidsitet,lsidsitet)
      implicit double precision (a-h,o-z)
      character str*100,sidsitet*25
      integer lsidsitet
      dimension pl(2)
      
      sidsitet='                         '

        pl(1)=0.0
        pl(2)=0.0
        do i=1,50
          if (str(i:i).eq.' ') then
             sidsitet=str(1:i-1)
             lsidsitet=i-1
             goto 81
          end if 
        end do
      

C   Latitude/longitude
 81   k=i  
      do il=1,2
      mpg1=0
      mp=0
      do m=k,k+30
      
      if ((str(m:m).eq.' ').and.(str(m+1:m+1).ne.' ')) then
                nsign=1
                m1=m+1
          
          if(str(m+1:m+1).eq.'-') then
          nsign=-1
          m1=m+2
            end if
            
          if(str(m+1:m+1).eq.'+') then
            nsign=+1
            m1=m+2
                end if
          
          do j=m1,m1+30
          if(str(j:j).eq.'.') then
          mp=j
          end if
          if(str(j:j).eq.':') then
          if(mpg1.eq.0) then
          mpg1=j
          else
          mpg2=j
          end if 
          end if
          
          if((str(j:j).eq.' ').and.(mp.eq.0).and.(mpg1.eq.0)) then
              idec=1
                mp=j
              mf=mp 
              goto 82
            end if
          if (mpg1.ne.0) then
              idec=0
           else
           idec=1
           end if
           
           if(str(j:j).eq.' ') then
           mf=j-1
           goto 82
           end if
           
           end do
           end if
           
           end do
   82 if (idec.eq.1) then
         nd=mf-mp
         ne=mp-m1
         p1=0
         
         do i=1,ne
            p1=p1+(ichar(str(m1+i-1:m1+i-1))-48)*10**(ne-i)
            end do
            
            p12=0
            do i=1,nd
            p12=p12+(ichar(str(mp+i:mp+i))-48)*(10**(nd-i))
            end do
            p12=p12/10.0**(nd)
            pl(il)=(p1+p12)*nsign
            k=mf+1
            else
            ndeg=mpg1-m1
            nmin=mpg2-mpg1-1
            if(mp.ne.0) then
            nsecd=mf-mp
            nsec=mp-mpg2-1
            else
            nsecd=0
            nsec=mf-mpg2
            end if
            
            p11=0
            do i=1,ndeg
            p11=p11+(ichar(str(m1+i-1:m1+i-1))-48)*10**(ndeg-i)
            
            end do
            p12=0
            do i=1,nmin
            p12=p12+(ichar(str(mpg1+i:mpg1+i))-48)*(10**(nmin-i))
            end do
            p13=0
            do i=1,nsec
            p13=p13+(ichar(str(mpg2+i:mpg2+i))-48)*(10**(nsec-i))
            end do
            p14=0
            do i=1,nsecd
            p14=p14+(ichar(str(mp+i:mp+i))-48)*(10**(nsecd-i))
            end do
            pl(il)=(p11+p12/60+(p13+p14/10**nsecd)/3600)*nsign
            k=mf+1
            end if
            end do
            p1=pl(1)
            p2=pl(2)
            return
            end

C   ***************************************************************************************
C     Subroutine Outputheaders.
C     Creates output files and writes headers
      
      SUBROUTINE Outputheaders(out,nout,th,statf,
     +           errov)
      implicit double precision (a-h,o-z)
      character file3*30, file4*30,file7*30,file8*30,file10*30
      character out*25,statf*1,overw*1
      integer nout
      
      write(file3,301) out(1:nout),'_INV_',th,'.dat'
  301        format(a,a5,f4.2,a4)

        open(unit=3,file=file3,status='new',err=303)
          goto 304
  303 write(6,305) out(1:nout)
  305 format('One or more output files (root ',a,') might already 
     +exist',/,' overwrite them (Y/N)? [Default: Y]')
        read(5,'(a)') overw
          if (overw.eq.'') then
          overw='Y'
          end if

        if((overw.eq.'n').or.(overw.eq.'N')) then

        write(6,*) 'You chose NOT to overwrite existing files'
          write(6,*) 'Program will close'
          write(6,*) 'Change output root name in param.cfg'
        errov=1
          close(3)
          goto 1000
        end if

  304   close(3)

        open(unit=3,file=file3,status='unknown')

      write(3,302) (in,in=1,7),(in,in=1,7)
  302 format('Site',27x,'Longitude',8x,
     +'Latitude',9x,'F(Hz)',4x,'Per(s)',8x,2('I',i1,'(km/s)',3x),6x,
     +5('I',i1,9x),'Q',9x,7('errI',i1,8x),1x,'errQ')
     
      write(file4,41) out(1:nout),'_DIM_',th,'.dat'
  41  format(a,a5,f4.2,a4)
  
      open(unit=4,file=file4,status='unknown')
      write(4,42) (il,il,il=3,7)
   42 format('Site',27x,'Longitude',8x,
     +'Latitude',9x,'F(Hz)',4x,'Per(s)',7x,'DIM',7x,
     +5('st',i1,4x,'errst',i1,7x))

      write(file7,701) out(1:nout),'_ERR_',th,'.dat'
 701  format(a,a5,f4.2,a4)

      open(unit=7,file=file7,status='unknown')
      write(7,702) (kl,kl=1,8),(kl,kl=1,9),(kl,kl=1,8),(kl,kl=1,9)
  702 format('Site',27x,'Longitude',8x,
     +'Latitude',9x,'F(Hz)',4x,'log(F)',6x,'Per(s)',6x,'ROT',7x,
     +'DIM',2x,'I',i1,'(km/s)',5x,'I',i1,'(km/s)',7x,6('I',i1,10x),
     +'st',i1,8(9x,'st',i1),1x,8(6x,'errI',I1),9(5x,'errst',i1))
      
      if ((statf.eq.'y').or.(statf.eq.'Y')) then
      
      write(file8,801) out(1:nout),'_STATS_',th,'.dat'
  801        format(a,a7,f4.2,a4)
      open(unit=8,file=file8,status='unknown')
      write(8,802) (ni,ni,ni,ni,ni,ni=1,8),(ni,ni,ni,ni,ni,ni=3,9)
  802 format('Site',27x,'Longitude',8x,
     +'Latitude',9x,'F(Hz)',4x,'Per(s)'
     +,4x,8(6x,'I',i1,'True',7x,'I',i1,'err',7x,'I',i1,'Sta',7x,
     +'I',i1,'dev',6x,'I',i1,'bias'),5x,7
     +('St',i1,'True',6x,'St',i1,'err',6x,'St',i1,'Sta',6x,'St',i1,'dev'
     +,4x,'St',i1,'nrmbi',5x))
      lc=0
      else
      file8=' '
      end if

      write(file10,1001) out(1:nout),'_other_',th,'.dat'
 1001        format(a,a7,f4.2,a4)
      open(unit=10,file=file10,status='unknown')
      write(10,1002)
 1002 format('Site',27x,'Per(s)',7x,'DIM',1x,'STRIKE(º)',14x,
     +'Impxy_rot',19x,'Impyx_rot',9x,'twist(º)',3x,
     +'shear(º)',9x,'I7',8x,'skew',3x,'ph_s_skew')

 1000 return ! end of Outputheaders
      end

C   ***************************************************************************************
C     Subroutine Header14.
C     Creates header of file 14.

      SUBROUTINE Outputheader14(out,nout,th)

      implicit double precision (a-h,o-z)

      character file14*30
      character out*25
      integer nout

      write(file14,1401) out(1:nout),'_BANDCLASS_',th,'.dat'
 1401 format(a,a11,f4.2,a4)


      open(unit=14,file=file14,status='unknown')
      write(14,140)
 140  format('Site',27x,'Longitude',8x,
     +'Latitude',3x,'BAND',7x,'Tmin',8x,'Tmax',1x,'nper',2x,'DIM',2x,
     +'strike',3x,'errstr',6x,'twi',3x,'errtwi',4x,'shear',4x,'errsh',
     +2x,'cont',2x,'scale',2x,'strikecomp')   
 
      return

        end

C   ***************************************************************************************
C     Subroutine Header15
C     Writes the header of file 15 (one file for each band).

      Subroutine Header15(out,nout,nc,th,permin,permax)
      
      implicit double precision (a-h,o-z)
      character file15*30,out*25
      integer nc,laux,nout
      character aux*2

      if (nc.lt.10) then
      write(aux,11) nc
        laux=1
        else
      write(aux,12) nc
      laux=2
        endif

   11 format(i1)
   12 format(i2)


      file15='                              '
      write(file15,1501) out(1:nout),'_BAND',aux(1:laux),'_',th,'.dat'
 1501 format(a,a5,a2,a,f4.2,a4)
      open(unit=15,file=file15,status='unknown')

      write(15,*)
      write(15,1502) nc, permin, permax
 1502 format('Dimensionality results average for band ',i3,':',
     +g10.3,' s - ',g10.3,' s')
      write(15,*)

      write(15,150)
  150 format('Site',27x,'Longitude',8x,
     +'Latitude',3x,'nper',2x,'DIM',2x,'strike',3x,'errstr',6x,'twi',3x,
     +'errtwi',4x,'shear',4x,'errsh',
     +2x,'cont',1x,'determ',2x,'strikecomp')
     
      return
      end
      
      
C   ***************************************************************************************
      Subroutine WriteFile3(mp,k,plongi,plati,sidsite,f,Y,devY)
      implicit double precision (a-h,o-z)
      Dimension Y(10,mp,2),devY(10,mp,2)
      Character sidsite*25

        
      write(3,73) sidsite,plongi,plati,
     +f,1/f,(abs(Y(j,k,1)),j=1,8),
     +(devY(j,k,1),j=1,8)
   73 format(a,1x,2(g14.6,2x),2x,
     +2(g10.3,2x),2x,
     +2(e10.3,1x),6(f10.5,1x),8(e12.5,1x))
     
      return ! end of WriteFile3
      end
C   ***************************************************************************************
      Subroutine WriteFile4(mp,k,sidsite,plongi,plati
     +                     ,f,icas,yst,devst)
      implicit double precision (a-h,o-z)
      Dimension yst(10,mp,2),devst(10,mp,2)
      Character sidsite*25

      write(4,403)sidsite,plongi,plati,f,1/f,icas,
     +(yst(la,k,1),devst(la,k,2),la=3,7)
  403 format(a,1x,2(g14.6,2x),2x,
     +2(g10.3,2x),5x,i1,1x,10(f9.4,1x))
      return ! end of WriteFile4
      end
C   ***************************************************************************************
      
      Subroutine WriteFile7(mp,k,li,sidsite,plongi,
     +                     plati,f,rot,icas,y,devy,yst,devst)
      implicit double precision (a-h,o-z)
      dimension y(10,mp,2),devy(10,mp,2),yst(10,mp,2),devst(10,mp,2)
      Character sidsite*25

      write(7,703) sidsite,plongi,plati,
     +f,log10(f),1/F,rot,
     + icas,(y(i,k,li),i=1,8),((yst(i,k,li)),i=1,9),
     + (devy(i,k,li),i=1,8),((devst(i,k,2)),i=1,9)
  703 format(a,1x,2(g14.6,2x),2x,
     +4(g10.3,2x),i1,2x,2(g11.4,2x),
     +6(g10.3,2x),g11.4,2x,8(g10.3,2x),17(g10.3,1x))
      
      return 
      end

C   ***************************************************************************************
      Subroutine WriteFile8(mp,k,sidsite,plongi,plati,
     +                      f,y,devy,yst,devst,bias,biasst)
      implicit double precision (a-h,o-z)
      dimension y(10,mp,2),devy(10,mp,2),yst(10,mp,2),devst(10,mp,2)
      dimension bias(10),biasst(10)
      character sidsite*25
      
      write(8,802) sidsite,plongi,plati,
     +f,1/F,
     +(y(i,k,1),devY(i,k,1),y(i,k,2),devY(i,k,2),bias(i),i=1,8),
     +(yst(i,k,1),devst(i,k,1),yst(i,k,2),devst(i,k,2),biasst(i),
     +i=3,9)
  802 format(a,1x,2(g14.6,2x),2x,
     +2(g10.3,2x),1x,75(e11.4,1x))
      return 
      
      end

C   ***************************************************************************************
      Subroutine WriteFile14(mp,sidsite,plongi,plati,kl,per,
     +           ipin,i1,i2,maxi,strike,errstr,twi,errtwi,shear,errsh,
     +                       nst,scale)

      implicit double precision (a-h,o-z)
      dimension per(mp)
      character sidsite*25
      integer i1,i2,ipin
      
      write(14,142) sidsite,plongi,plati,kl,per(i1),per(i2),
     +              ipin,maxi,strike,errstr,twi,errtwi,shear,errsh,nst,
     +              scale,strike*(-1)
  142 format(a,1x,2(g14.6,2x),
     +1x,i3,2(f12.4),i4,1x,i3,1x,6(f9.4),2x,i2,2(f9.4))
     
      return
     
      end

C   ***************************************************************************************
      Subroutine Header9(sidsite,lsidsite,th,file9)

      implicit double precision (a-h,o-z)

      Character sidsite*25
      Character file9*30
      Integer lsidsite
      write(file9,91) sidsite(1:lsidsite),'_RES_',th,'.dat'
   91        format(a,a5,f4.2,a4)

      open(unit=9,file=file9,status='unknown')


      WRITE(9,*)
      write(9,92) sidsite(1:lsidsite)
   92 format('Writing main dimensionality analysis results for site ',
     +a)
      WRITE(9,*)
      write(9,*) 'Table with results for each frequency'
      WRITE(9,*)
      write(9,93) 
   93 format('Site',27x,'Longitude',8x,
     +'Latitude',9x,'F(Hz)',4x,'Per(s)',10x,'DIM',3x,'Strike(º)',
     +1x,'errStrike(º)',
     +5x,'Twist(º)',2x,'errTwist(º)',5x,'Shear(º)',
     +2x,'errShear(º)')

      return 
      end

C   ***************************************************************************************
      Subroutine WriteFile9(mp,k,sidsite,plongi,plati
     +                     ,f,icas,yst,devst,namedim)

      implicit double precision (a-h,o-z)
 
      Dimension yst(10,mp,2),devst(10,mp,2)
      Character sidsite*25,namedim(8)*9
     
C  Dimensionality type into character string:
C  argument = icas + 1.

      namedim(1)='UNDETERM.'
      namedim(2)='       1D'
      namedim(3)='       2D'
      namedim(4)='  3D/2Dtw'
      namedim(5)='    3D/2D'
      namedim(6)='       3D'
      namedim(7)='3D/2Ddiag'
      namedim(8)='  3D/1D2D'    

      if ((icas.lt.2).or.(icas.gt.4)) then

         write(9,901)sidsite,plongi,plati,f,1/f,namedim(icas+1)
  901 format(a,1x,2(g14.6,2x),2x,2(e10.3,2x),a9)
     
      end if
  

      if (icas.eq.2)  then
      
         st=((yst(3,k,1)+yst(4,k,1)))/2
         errst=((devst(3,k,2)+devst(4,k,2)))/2

         write(9,902)sidsite,plongi,plati,f,1/f,namedim(icas+1),
     +st,errst
  
  902    format(a,1x,2(g14.6,2x),2x,2(e10.3,2x),a9,3x,
     +2(f9.4,4x))
     
      end if
      
      if ((icas.eq.3).or.(icas.eq.4))  then
      
         st=yst(5,k,1)
         errst=devst(5,k,2)

         ti=yst(8,k,1)
         eti=devst(8,k,2)

         se=yst(9,k,1)
         ese=devst(9,k,2)

         write(9,903)sidsite,plongi,plati,f,1/f,namedim(icas+1),
     +st,errst,ti,eti,se,ese
  903    format(a,1x,2(g14.6,2x),2x,2(e10.3,2x),a9,3x,
     +6(f9.4,4x))
      
      end if
 
      return ! end of WriteSiteRes
      end

C   ***************************************************************************************
      Subroutine Header9bis

      implicit double precision (a-h,o-z)


      write(9,*)
      write(9,'(a)') 'Dimensionality results averaged in bands'
      write(9,*)

      write(9,190)
  190 format('Site',27x,'Longitude',8x,
     +'Latitude',3x,'BAND',7x,'Tmin',8x,'Tmax',1x,'nper',7x,'DIM',4x,
     +'strike',3x,'errstr',6x,'twi',3x,'errtwi',4x,'shear',4x,'errsh',
     +2x,'cont',2x,'scale',2x,'strikecomp')  


      return

      end

C   ***************************************************************************************
      Subroutine WriteFile9bis(mp,sidsite,plongi,plati,kl,per,
     +           ipin,i1,i2,namedim,maxi,strike,errstr,twi,errtwi,
     +           shear,errsh,nst,scale)

      implicit double precision (a-h,o-z)
      dimension per(mp)
      character sidsite*25,namedim(8)*9
      integer i1,i2,ipin

      write(9,192) sidsite,plongi,plati,kl,per(i1),per(i2),ipin,
     +namedim(maxi+1),strike,errstr,twi,errtwi,
     +shear,errsh,nst,scale,strike*(-1)
  192 format(a,1x,2(g14.6,2x),
     +1x,i3,2(f12.4),i3,3x,a9,1x,6(f9.4),2x,i2,2(f9.4))

        return 
        end

C   ***************************************************************************************
      Subroutine CompWriteFile10(mp,sidsite,k,f,
     +                          zr,zi,fi,pnu,y,yst,icas,namedim)

      implicit double precision (a-h,o-z)
      dimension zr(2,2,mp),zi(2,2,mp),fi(4,mp),pnu(4,mp)
      dimension y(10,mp,2),yst(10,mp,2)
      dimension rm(2,2)
      character sidsite*25,sim1*3,sim2*3,namedim(8)*9
      integer icas

C     Strike determined by the invariants:      
      strike=0.0
      
      if (icas.eq.2) then   
         strike=((yst(3,k,1)+yst(4,k,1)))/2
      end if

      if ((icas.eq.3).or.(icas.eq.4))  then     
         strike=yst(5,k,1)
      endif

C     rotation matrix

      str=strike*acos(-1.)/180.   

      rm(1,1)=cos(str)
      rm(1,2)=sin(str)
      rm(2,1)=-sin(str)
      rm(2,2)=cos(str)

C     Rotated impedances:

C     Non-diag:

      zxyre_rot=rm(1,1)*rm(2,1)*zr(1,1,k)+
     +rm(1,1)*rm(2,2)*zr(1,2,k)+
     +rm(1,2)*rm(2,1)*zr(2,1,k)+
     +rm(1,2)*rm(2,2)*zr(2,2,k)

      zxyim_rot=rm(1,1)*rm(2,1)*zi(1,1,k)+
     +rm(1,1)*rm(2,2)*zi(1,2,k)+
     +rm(1,2)*rm(2,1)*zi(2,1,k)+
     +rm(1,2)*rm(2,2)*zi(2,2,k)

      if (zxyim_rot.ge.0) then
         sim1=' + '
         else
         sim1=' - '
      endif

      zyxre_rot=rm(1,1)*rm(2,1)*zr(1,1,k)+
     +rm(1,2)*rm(2,1)*zr(1,2,k)+
     +rm(1,1)*rm(2,2)*zr(2,1,k)+
     +rm(1,2)*rm(2,2)*zr(2,2,k)

      zyxim_rot=rm(1,1)*rm(2,1)*zi(1,1,k)+
     +rm(1,2)*rm(2,1)*zi(1,2,k)+
     +rm(1,1)*rm(2,2)*zi(2,1,k)+
     +rm(1,2)*rm(2,2)*zi(2,2,k)

      if (zyxim_rot.ge.0) then
         sim2=' + '
         else
         sim2=' - '
      endif



C     Bahr parameters: Skew (b1) and phase sensitive skew (b3)

      b1=sqrt(fi(1,k)**2+pnu(1,k)**2)/sqrt(fi(4,k)**2+pnu(4,k)**2)
 
      b3=sqrt(abs(fi(3,k)*pnu(2,k)-fi(2,k)*pnu(3,k)-
     +fi(1,k)*pnu(4,k)+fi(4,k)*pnu(1,k)))/Yinv12(fi(4,k),pnu(4,k))


      if ((icas.le.1).or.(icas.ge.6)) then

         write(10,1011) sidsite,1/f,namedim(icas+1),
     +   b1,b3

 1011    format(a,1x,e10.3,2x,a9,2x,'    ---',3x,
     +'     ---------------------',2x,
     +'     ---------------------',2x,
     +'    ------',1x,'    -----',3x,
     +' --------',2x,
     +g10.3,2x,g10.3)

      endif

      if (icas.eq.2) then

         write(10,1012) sidsite,1/f,namedim(icas+1),strike,
     +   '(',zxyre_rot,sim1,abs(zxyim_rot),' i)',
     +   '(',zyxre_rot,sim2,abs(zyxim_rot),' i)',
     +   b1,b3

 1012    format(a,1x,e10.3,2x,a9,2x,f7.2,3x,
     +a1,f10.4,a3,f9.4,a3,2x,
     +a1,f10.4,a3,f9.4,a3,2x,
     +'    ------',1x,'    -----',3x,
     +' --------',2x,
     +g10.3,2x,g10.3)

      endif


      if ((icas.eq.3).or.(icas.eq.4)) then

         write(10,1013) sidsite,1/f,namedim(icas+1),strike,
     +   yst(8,k,1),yst(9,k,1),
     +   b1,b3

 1013    format(a,1x,e10.3,2x,a9,2x,f7.2,3x,
     +'     ---------------------',2x,
     +'     ---------------------',2x,
     +1x,f9.4,1x,f9.4,3x,
     +' --------',2x,
     +g10.3,2x,g10.3)

      endif

      if (icas.eq.5) then

         write(10,1015) sidsite,1/f,namedim(icas+1),
     +   y(7,k,1),b1,b3

 1015    format(a,1x,e10.3,2x,a9,2x,'    ---',3x,
     +'     ---------------------',2x,
     +'     ---------------------',2x,
     +'    ------',1x,'    -----',3x,
     +g10.3,1x,g10.3,2x,g10.3)
      
      endif

      return
      
      end