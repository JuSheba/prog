gf=gfortran -O3
files=$(wildcard *.f90)
prog: $(patsubst %.f90, %.o, $(files))
	$(gf) $^ -o $@
solve_stiff_DE.o: param_func_mod.o stiff_DE_mod.o
%.o : %.f90
	$(gf) -c $<
stiff_DE_mod.o: elimination_mod.o param_func_mod.o methods_adams_mod.o
%.o : %.f90
	$(gf) -c $<
methods_adams_mod.o: param_func_mod.o rk_mod.o coeff_adams_mod.o newton_mod.o
%.o : %.f90
	$(gf) -c $<
newton_mod.o: elimination_mod.o
%.o : %.f90
	$(gf) -c $<
rk_mod.o: param_func_mod.o
%.o : %.f90
	$(gf) -c $<
exe: prog
	./$< rk
	./$< ae
	./$< ai
clean:
	rm -f *.o prog *.mod
.PHONY: clean, exe
