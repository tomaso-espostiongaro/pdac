include Machine

all: comm.a pdac.x

FOBJS= atmosphere.o bdry.o control.o decomp.o dens.o dimensions.o drag.o \
environment.o eosg.o eosl.o ftem.o flux_u.o flux_v.o flux_w.o \
flux_sc.o flux_m.o gas.o grid.o hcapgs.o hotc.o hvs.o htilde.o \
kinds.o immb.o indijk.o input.o interp.o io.o iter.o kb07ad.o limiters.o\
letter.o matrix.o nondim.o outp.o parallel.o particles.o \
pdac.o press.o prog.o reactions.o roughness.o setup.o \
subscr.o temp.o tilde.o time.o turbo.o types.o \
velocity.o visc.o ygas.o iotk_module.o \
$(SYSOBJ) 

PPFOBJS= atmosphere.o bdry.o control.o decomp.o dens.o dimensions.o drag.o \
environment.o eosg.o eosl.o ftem.o flux_u.o flux_v.o flux_w.o \
flux_sc.o flux_m.o gas.o grid.o hcapgs.o hotc.o hvs.o htilde.o \
kinds.o immb.o indijk.o input.o io.o iter.o kb07ad.o limiters.o\
letter.o matrix.o nondim.o outp.o parallel.o particles.o \
press.o prog.o reactions.o roughness.o setup.o \
subscr.o temp.o tilde.o time.o turbo.o types.o \
velocity.o visc.o ygas.o filter.o pp.o iotk_module.o \
$(SYSOBJ) 

pdac.x: $(FOBJS) comm/comm.a
	$(LINKER) -o pdac.x $(MPFFLAGS) $(LINKFLAGS) $(FOBJS) comm/comm.a $(LIBS)

pp.x: $(PPFOBJS) comm/comm.a
	$(LINKER) -o pp.x $(MPFFLAGS) $(LINKFLAGS) $(PPFOBJS) comm/comm.a $(LIBS)

clean:
	rm -f *.o *.mod *.a core core.* *.stb
	(cd comm; make clean)

comm/comm.a comm.a:
	(cd comm; make comm.a)

util :
	(cd utility; cc -O2 -o moduledep.x moduledep.c)

include .dependencies
