BIN = mc_s mc_p ge_s ge_p ge_pg
MODS=readMatrix.o

MODSRC=$(patsubst %.o,%.f,$(MODS))
MODMPI=$(patsubst %.o,%mpi.o,$(MODS))

DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

mc_s: mc_s.o $(MODS)
	gfortran -o $@ $^ $(LFLAGS)
mc_p: mc_p.o $(MODMPI)
	mpifort -o $@ $^ $(LFLAGS)
ge_s: ge_s.o $(MODS)
	gfortran -o $@ $^ $(LFLAGS)
ge_p: ge_p.o $(MODMPI)
	mpifort -o $@ $^ $(LFLAGS)
ge_pg: ge_pg.o $(MODMPI)
	mpifort -o $@ $^ $(LFLAGS)

mc_s.o: mc_s.f95 $(MODS)
	gfortran -o $@ $< $(CFLAGS) -O3
ge_s.o: ge_s.f95 $(MODS)
	gfortran -o $@ $< $(CFLAGS) -O3
mc_p.o: mc_p.f95 $(patsubst %.o,%mpi.o,$(MODS))
	mpifort -o $@ $< $(CFLAGS) -O3
ge_p.o: ge_p.f95 $(patsubst %.o,%mpi.o,$(MODS))
	mpifort -o $@ $< $(CFLAGS) -O3
ge_pg.o: ge_pg.f95 $(patsubst %.o,%mpi.o,$(MODS))
	mpifort -o $@ $< $(CFLAGS) -O3

$(MODS) : $(MODSRC)
	gfortran -o $@ $< $(CFLAGS)

$(MODMPI) : $(MODSRC)
	mpifort -o $@ $< $(CFLAGS)

.PHONY: clean
clean:
	rm -f *.o *~ core $(BIN) *.mod
