include Machine

FOBJS= atmosphere.o bdry.o coll.o dens.o dimensions.o distr.o\
 drag.o environment.o eosg.o eosl.o ftem.o fluxes.o gas.o grid.o \
hcapgs.o hotc.o hvs.o htilde.o io.o iter.o kb07ad.o\
letter.o matrix.o outp.o parallel.o particles.o\
pdac2d.o press.o prog.o reactions.o roughness.o setc.o setup.o \
subscl.o temp.o tilde.o time.o turbo.o types.o\
velocity.o visc.o ygas.o $(SYSOBJ)

include .dependencies

pdac2d.x_sc: $(FOBJS) comm/comm.a
	$(LINKER) -o pdac2d.x_sc $(LINKFLAGS) $(FOBJS) comm/comm.a $(LIBS)

clean:
	rm -f *.o *.mod *.a core core.*
	(cd comm; make clean)

comm/comm.a comm.a:
	(cd comm; make comm.a)

include Rules
