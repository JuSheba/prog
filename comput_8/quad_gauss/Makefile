gf=gfortran -O3
files=$(wildcard *.f90)
prog: $(patsubst %.f90, %.o, $(files))
	$(gf) $^ -o $@
main.o: quad_gauss.o
%.o : %.f90
	$(gf) -c $<
quad_gauss.o: legendre_mod.o elimination_mod.o
%.o : %.f90
	$(gf) -c $<
legendre_mod.o: bernoulli_mod.o
%.o : %.f90
	$(gf) -c $<
bernoulli_mod.o: gorner_mod.o
%.o : %.f90
	$(gf) -c $<
exe: prog
	./$<
clean:
	rm -f *.o prog *.mod
.PHONY: clean, exe
