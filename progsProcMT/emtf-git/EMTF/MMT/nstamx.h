      integer nstamx,nchmx,ntmx,nstamx2,nchemx
      parameter(nstamx = 8,nchmx = 12,ntmx=nstamx*nchmx)
      parameter(nstamx2 = 2*nstamx,nbmax=50)
ccc  set nchemx to largest n umber of predicted channels in
ccc    any TF group (by default each TF group is a station ...
ccc    but this can be changed (e.g., all channels in one group)
ccc    in this case nchemx needs to be at least nt - 2 where nt is 
ccc    the total number of channels
       parameter(nchemx=12,ngrpmx = 8)
