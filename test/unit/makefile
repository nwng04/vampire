#===================================================================
#
#             Makefile for Vampire Test Suite
#
#===================================================================

# Compilers
GCC=g++

# LIBS
LIBS=-lstdc++

# Flags
GCC_CFLAGS=-O3 -std=c++17 -I../../hdr/ -I../../src/ -DDEBUG

# Objects
TEST_OBJECTS= \
obj/unit_tests.o \
obj/utility/units_test.o \
obj/utility/utility_test.o\
obj/utility/spin_temperature_test.o \
obj/spinlattice/interface_test.o \
obj/spinlattice/suzuki_trotter_test.o \
obj/spinlattice/suzuki_trotter_stubs.o


VAMPIRE_OBJECTS= \
../../obj/main/githash.o \
../../obj/main/version.o \
../../obj/main/material.o \
../../obj/utility/errors.o \
../../obj/utility/units.o \
../../obj/vio/data.o \
../../obj/vio/globalio.o \
../../obj/vio/timestamp.o\
../../obj/vio/check.o\
../../obj/vio/string_to_x.o\
../../obj/constants/constants.o\
../../obj/spinlattice/temperatures.o\
../../obj/spinlattice/suzuki-trotter.o\
../../obj/spinlattice/data.o \
../../obj/data/atoms.o\
../../obj/spinlattice/interface.o


EXECUTABLE=unit_tests

all: $(TEST_OBJECTS) $(VAMPIRE_OBJECTS) gcc

# Serial Targets
gcc: $(TEST_OBJECTS) $(VAMPIRE_OBJECTS)
	$(GCC) $(TEST_OBJECTS) $(VAMPIRE_OBJECTS) $(GCC_CFLAGS) $(LIBS) -o $(EXECUTABLE)

$(TEST_OBJECTS): obj/%.o: src/%.cpp
	$(GCC) -c -o $@ $(GCC_CFLAGS) $<

# compile vampire objects to vampire object folder
$(VAMPIRE_OBJECTS): ../../obj/%.o: ../../src/%.cpp
	$(GCC) -c -o $@ $(GCC_CFLAGS) $<

clean:
	@rm -f obj/*.o
	@rm -f ../../obj/*.o

purge:
	@rm -f obj/*.o
	@rm -f $(EXECUTABLE)
