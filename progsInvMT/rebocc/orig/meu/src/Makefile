
OBJ = sinv2d.o getdata.o outdata.o modelsetup.o \
      mcov.o setupop.o sens2d.o representer.o\
      fwd2d.o fwd1d.o search.o\
      utilfunct.o sort_median.o

rebocc : $(OBJ)
	gfortran -g $(OBJ) -o rebocc -m64 -lacml_mp

sinv2d.o : sinv2d.f
	gfortran -g -O2 -c sinv2d.f

getdata.o : getdata.f
	gfortran -g -O2 -c getdata.f

outdata.o : outdata.f
	gfortran -g -O2 -c outdata.f

modelsetup.o : modelsetup.f
	gfortran -g -O2 -c modelsetup.f

mcov.o : mcov.f
	gfortran -g -O2 -c mcov.f

setupop.o : setupop.f
	gfortran -g -O2 -c setupop.f

sens2d.o : sens2d.f
	gfortran -g -O2 -c sens2d.f

representer.o : representer.f
	gfortran -g -O2 -c representer.f

search.o : search.f
	gfortran -g -O2 -c search.f

fwd2d.o : fwd2d.f
	gfortran -g -O2 -c fwd2d.f

fwd1d.o : fwd1d.f
	gfortran -g -O2 -c fwd1d.f

utilfunct.o : utilfunct.f
	gfortran -g -O2 -c utilfunct.f

sort_median.o : sort_median.f
	gfortran -g -O2 -c sort_median.f
