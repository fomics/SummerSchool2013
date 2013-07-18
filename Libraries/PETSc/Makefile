#
# make all builds both the GPU and host versions versions
#
all:	FEsolver.exe
        
%.o: %.c
	cc -c -O3 $< -I.
#	cc -c -O3 $< -I. -I/opt/petsc/3.3.00/real/GNU/46/include -I/opt/petsc/3.3.00/real/GNU/46/interlagos/include

FEsolver.exe: FEsolver.o 
	cc -O3 -o FEsolver.exe FEsolver.o

FEsolver_ex.exe: FEsolver_ex.o 
	cc -O3 -o FEsolver_ex.exe FEsolver_ex.o

FEsolver_complete.exe: FEsolver_complete.o
	cc -O3 -o FEsolver_complete.exe FEsolver_complete.o

clean:
	rm -f *.o *.exe 
