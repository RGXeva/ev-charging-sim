# Build with Open MPI + OpenMP.
# macOS: Apple clang has no OpenMP, so mpicc wraps Homebrew GCC (override with OMPI_CC=...).
# Linux / SLURM clusters: use an Open MPI whose mpicc wraps GCC (module load openmpi ...).

UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
  export OMPI_CC ?= $(firstword $(notdir $(wildcard /opt/homebrew/bin/gcc-1[0-9])) gcc)
endif

MPICC    ?= mpicc
CFLAGS   ?= -O2 -std=c11 -Wall -Wextra -Wpedantic -fopenmp
LDLIBS   := -lm
SRC      := src/main.c src/util.c src/node.c src/base.c src/base_logic.c
HDR      := src/common.h
LAUNCHER ?= mpirun --oversubscribe

all: bin/wsn bin/test_candidates

bin/wsn: $(SRC) $(HDR) | bin
	$(MPICC) $(CFLAGS) -o $@ $(SRC) $(LDLIBS)

# MPI-free unit test of the base-station candidate search
bin/test_candidates: tests/test_candidates.c src/base_logic.c $(HDR) | bin
	$(MPICC) $(CFLAGS) -Isrc -o $@ tests/test_candidates.c src/base_logic.c $(LDLIBS)

# Debug build with AddressSanitizer + UBSan (Linux; launch with mpirun, not srun)
bin/wsn_asan: $(SRC) $(HDR) | bin
	$(MPICC) -O1 -g -std=c11 -Wall -Wextra -fopenmp -fsanitize=address,undefined -fno-omit-frame-pointer \
	    -o $@ $(SRC) $(LDLIBS)

asan: bin/wsn_asan

# Unit test + the full integration suite (every version pair, shapes, scenarios, bad input)
test: all
	./bin/test_candidates
	python3 -I tests/run_tests.py --launcher "$(LAUNCHER)"

# Unit test + a short subset of the integration suite
test-quick: all
	./bin/test_candidates
	python3 -I tests/run_tests.py --launcher "$(LAUNCHER)" --quick

# A 3 x 3 grid + base station for 20 s, printing every alert and the summary
demo: bin/wsn
	LAUNCHER="$(LAUNCHER)" bash scripts/demo_local.sh

bin:
	mkdir -p bin

clean:
	rm -rf bin out

.PHONY: all asan test test-quick demo clean
