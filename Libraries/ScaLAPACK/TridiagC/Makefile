#
# make all builds both the GPU and host versions versions
#
all:	pddttrdrv.exe
        
%.o: %.c
	cc -c -O3 $< -I. 

pddttrdrv.exe: pddttrdrv.o
	cc -O3 -o pddttrdrv.exe pddttrdrv.o

pddttrdrv_complete.exe: pddttrdrv_complete.o
	cc -O3 -o pddttrdrv_complete.exe pddttrdrv_complete.o

clean:
	rm -f *.o *.exe *.mod *.MOD *.out
