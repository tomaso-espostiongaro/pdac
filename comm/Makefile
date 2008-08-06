include ../Machine
MODULES=

COMM : barrier.o startup.o hangup.o bcast.o localdim.o gather.o scatter.o sendrecv.o \
	reduce.o error.o wtime.o myflush.o cptimer.o isend.o recv.o wait.o irecv.o \
        floatchk.o arraychk.o memstat.o commutil.o
#	ar $(ARFLAGS) $@ $?

clean : 
	rm -rf *.a *.o core*

