include ../Machine

comm.a : startup.o hangup.o bcast.o localdim.o gather.o scatter.o sendrecv.o \
	reduce.o error.o
	ar rv $@ $?

clean : 
	rm -rf *.a *.o core*

include ../Rules
