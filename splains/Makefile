test=tester.f90
files=$(filter-out $(test), $(wildcard *.f90))
obj=$(patsubst %.f90, %.o, $(files))

res: result

result: main
main: $(obj) $(test)
	gfortran $^  -o  $@
%.o: %.f90
	gfortran -c $<  -o $@
$(test):
	gfortran $@ -o $(test: .f90 = .o)
clean:
	rm -f *.o *.mod *~ main
