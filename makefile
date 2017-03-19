PNAME=equilibrium_QMC
CC=g++
CFLAGS=-I.
LIBS=
DEPS=MC_main.hpp MC_random.hpp MC_coord.hpp MC_finiteTemperatureEquilibriumQMC.hpp MC_histogram.hpp MC_heisenbergEquilibriumQMC.hpp MC_isingEquilibriumMC.hpp
OBJ=MC_main.o MC_random.o MC_finiteTemperatureEquilibriumQMC.o MC_histogram.o MC_heisenbergEquilibriumQMC.o MC_isingEquilibriumMC.o

all: $(PNAME)

%.o: %.cpp $(DEPS)
	$(CC) -O3 -std=c++0x -c -o $@ $< $(CFLAGS) $(LIBS)

$(PNAME): $(OBJ)
	$(CC) -O3 -std=c++0x -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm  *.o *~
