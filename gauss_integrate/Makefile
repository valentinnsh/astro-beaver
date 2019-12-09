test=tester.f90
files=$(filter-out $(test), $(wildcard *.f90))
obj=$(patsubst %.f90, %.o, $(files))

res: result

result: main
main: horner_and_bernoulli.o file_ops.o $(obj) $(test)
	gfortran $^  -o  $@
file_ops.o: file_ops.f90
	gfortran -c $< -o $@
horner_and_bernoulli.o: horner_and_bernoulli.f90
	gfortran -c $< -o $@
%.o: %.f90
	gfortran -c $<  -o $@
$(test):
	gfortran $@ -o $(test: .f90 = .o)
clean:
	rm -f *.o *.mod *~ main
