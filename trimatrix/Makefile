files=$(wildcard *.f90)
obj=$(patsubst %.f90, %.o, $(files))

res: result clean

result: main
	./main
main: $(obj)                                                          
	gfortran $^ -o $@
%.o: %.f90
	gfortran -c $< -o $@
clean:
	rm -f *.o main
