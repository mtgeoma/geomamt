
OBJ = sinv2d.o getdata.o outdata.o modelsetup.o \
      mcov.o setupop.o sens2d.o representer.o\
      fwd2d.o fwd1d.o search.o\
      utilfunct.o\
      LAPACK_R8.o LAPACK_C16.o BLAS1.o BLAS2_3.o 

FLAGS= -fomit-frame-pointer -mfpmath=387 -O2 -falign-loops=4 -mcmodel=medium -m64

rebocc : $(OBJ) 
	gfortran $(FLAGS) $(OBJ) -o rebocc 

sinv2d.o : sinv2d.f 
	gfortran $(FLAGS) -c sinv2d.f 

getdata.o : getdata.f 
	gfortran $(FLAGS) -c getdata.f 

outdata.o : outdata.f 
	gfortran $(FLAGS) -c outdata.f 

modelsetup.o : modelsetup.f 
	gfortran $(FLAGS) -c modelsetup.f 

mcov.o : mcov.f 
	gfortran $(FLAGS) -c mcov.f 

setupop.o : setupop.f 
	gfortran $(FLAGS) -c setupop.f 

sens2d.o : sens2d.f 
	gfortran $(FLAGS) -c sens2d.f 

representer.o : representer.f 
	gfortran $(FLAGS) -c representer.f 

search.o : search.f 
	gfortran $(FLAGS) -c search.f 

fwd2d.o : fwd2d.f 
	gfortran $(FLAGS) -c fwd2d.f 

fwd1d.o : fwd1d.f 
	gfortran $(FLAGS) -c fwd1d.f 

utilfunct.o : utilfunct.f 
	gfortran $(FLAGS) -c utilfunct.f 

BLAS1.o : BLAS1.f
	gfortran $(FLAGS) -c BLAS1.f
 
BLAS2_3.o : BLAS2_3.f
	gfortran $(FLAGS) -c BLAS2_3.f
 
LAPACK_R8.o : LAPACK_R8.f
	gfortran $(FLAGS) -c LAPACK_R8.f
 
LAPACK_C16.o : LAPACK_C16.f
	gfortran $(FLAGS) -c LAPACK_C16.f
