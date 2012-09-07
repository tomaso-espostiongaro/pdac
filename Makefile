include Machine

FOBJS= arrfil.o atmosphere.o blunt.o bdry.o control.o decomp.o dens.o \
dimensions.o dome.o drag.o environment.o eosg.o eosl.o ftem.o \
flux_u.o flux_v.o flux_w.o flux_sc.o flux_m.o gas.o ghost.o grid.o \
hcapgs.o hotc.o hvs.o htilde.o kinds.o immb.o indijk.o input.o \
interp.o inoutflow.o io_parallel.o io_restart.o iter.o kb07ad.o \
limiters.o letter.o matrix.o\
outp.o parallel.o particles.o press.o prog.o reactions.o rannum.o\
residuals.o roughness.o samplert.o setup.o sink.o subscr.o temp.o tilde.o time.o \
topo.o \
turbo.o types.o velocity.o vent.o visc.o ygas.o iotk_module.o io_files.o

PPFOBJS = derived.o filter.o mean_fields.o masspart.o massflux.o massgsedim.o \
process.o postin.o postout.o postp.o postvar.o sampling.o sections.o

all: COMM pdac.x postp.x

pdac.x: $(FOBJS) pdac.o COMM 
	$(LINKER) -o pdac.x $(MPFFLAGS) $(LINKFLAGS) $(FOBJS) pdac.o comm/*.o $(LIBS)

postp.x: $(FOBJS) $(PPFOBJS) COMM
	$(LINKER) -o postp.x $(MPFFLAGS) $(LINKFLAGS) $(FOBJS) $(PPFOBJS) comm/*.o $(LIBS)

tstcomm.x: testcommlib.o COMM
	$(LINKER) -o tstcomm.x $(MPFFLAGS) $(LINKFLAGS) testcommlib.o comm/*.o $(LIBS)

clean:
	rm -f *.o *.mod *.a core core.* *.stb
	(cd comm; make clean)

COMM:
	(cd comm; make )

util :
	(cd utility; cc -O2 -o moduledep.x moduledep.c)

include .dependencies
