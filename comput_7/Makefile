gf=gfortran -O3
files=$(wildcard *.f90)
prog: $(patsubst %.f90, %.o, $(files))
	$(gf) $^ -o $@
fourier.o: fourier_module.o
%.o : %.f90
	$(gf) -c $<
exe: prog
	./$<
clean:
	rm -f *.o prog *.mod
.PHONY: clean, exe
