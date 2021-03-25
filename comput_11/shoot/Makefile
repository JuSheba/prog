gf=gfortran -O3
files=$(wildcard *.f90)
prog: $(patsubst %.f90, %.o, $(files))
	$(gf) $^ -o $@
boundary_value_prblm.o: shooting_method_mod.o param_func.o
%.o : %.f90
	$(gf) -c $<
shooting_method_mod.o:  newton_mod.o param_func.o rk_mod.o
%.o : %.f90
	$(gf) -c $<
newton_mod.o: elimination_mod.o
%.o : %.f90
	$(gf) -c $<
rk_mod.o: param_func.o
%.o : %.f90
	$(gf) -c $<
exe: prog
	./prog
clean:
	rm -f *.o prog *.mod
.PHONY: clean, exe
