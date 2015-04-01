
OBJ = sinv2d.o getdata.o outdata.o modelsetup.o \
      mcov.o setupop.o sens2d.o representer.o\
      fwd2d.o fwd1d.o search.o\
      utilfunct.o\
      LAPACK_R8.o LAPACK_C16.o BLAS1.o BLAS2_3.o 

rebocc : $(OBJ) 
	g77 -g $(OBJ) -o rebocc 

sinv2d.o : sinv2d.f 
	g77 -g -O -c sinv2d.f 

getdata.o : getdata.f 
	g77 -g -O -c getdata.f 

outdata.o : outdata.f 
	g77 -g -O -c outdata.f 

modelsetup.o : modelsetup.f 
	g77 -g -O -c modelsetup.f 

mcov.o : mcov.f 
	g77 -g -O -c mcov.f 

setupop.o : setupop.f 
	g77 -g -O -c setupop.f 

sens2d.o : sens2d.f 
	g77 -g -O -c sens2d.f 

representer.o : representer.f 
	g77 -g -O -c representer.f 

search.o : search.f 
	g77 -g -O -c search.f 

fwd2d.o : fwd2d.f 
	g77 -g -O -c fwd2d.f 

fwd1d.o : fwd1d.f 
	g77 -g -O -c fwd1d.f 

utilfunct.o : utilfunct.f 
	g77 -g -O -c utilfunct.f 

BLAS1.o : BLAS1.f
	g77 -g -O -c BLAS1.f
 
BLAS2_3.o : BLAS2_3.f
	g77 -g -O -c BLAS2_3.f
 
LAPACK_R8.o : LAPACK_R8.f
	g77 -g -O -c LAPACK_R8.f
 
LAPACK_C16.o : LAPACK_C16.f
	g77 -g -O -c LAPACK_C16.f
