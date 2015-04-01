c****** number of channels = nch
ccc     set this to the exact number of channels in the files
ccc     to be formatted
ccc     ALSO: change chid
      parameter (nch = 5,nblkmx=20000)
      parameter(nchp1 = nch+1, ntotmx = 1000)
      character*10 chid(nch)
cc    5 channel MT
      data chid/'Hx','Hy','Hz','Ex','Ey'/
