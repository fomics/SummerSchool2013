#
# make all builds both the GPU and host versions versions
#
all:	FEsolver.exe
        
%.o: %.c
	cc -c -O3 $< -I. 

FEsolver.exe: FEsolver.o 
	cc -O3 -o FEsolver.exe FEsolver.o

FEsolver_complete.exe: FEsolver_complete.o
	cc -O3 -o FEsolver_complete.exe FEsolver_complete.o

clean:
	rm -f *.o *.exe 
