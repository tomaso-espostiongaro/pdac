include ../Machine
MODULES=

comm.a : startup.o hangup.o bcast.o localdim.o gather.o scatter.o sendrecv.o \
	reduce.o error.o wtime.o myflush.o cptimer.o isend.o recv.o wait.o irecv.o
	ar rv $@ $?

clean : 
	rm -rf *.a *.o core*

include ../Rules
