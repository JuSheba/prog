  COMP   :=gfortran
  PATTERN:=*.f95
  SOURCE :=$(wildcard $(PATTERN))
  OBJECT :=$(patsubst %.f95, %.o, $(SOURCE))
    main : $(OBJECT)
	   $(COMP) $^ -o $@
 %.o %.mod : %.f95
	   $(COMP) -c -Wall $<
 FirstProg.o: my_prec.mod

   my_prec.o my_prec.mod : my_prec.f95
	 $(COMP) -c my_prec.f95
   clear :
	 rm -f *.o *.mod main res; touch res
	
  result: main input  
	  ./main < input > result

