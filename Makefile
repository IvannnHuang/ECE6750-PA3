CC ?= gcc
# CFLAGS ?=
CFLAGS ?= -Wall -Wextra -O3 -march=native -ffast-math -funroll-loops
LDFLAGS ?= -lm

IMPL ?= pthreads
MODE ?= parallel
matrix ?= small.mtx
p ?= 8
test ?= 0
PARTITION ?= ece6750

ifeq ($(test),1)
TEST_FLAG = -DTEST
else
TEST_FLAG =
endif

.PHONY: all full clean pspgemm serial

all: full

ifeq ($(MODE),single)
RUN_BIN = serial
RUN_P =
RUN_TAG = serial_$(basename $(notdir $(matrix)))
else
RUN_BIN = pspgemm
RUN_P = $(p)
RUN_TAG = $(IMPL)_$(basename $(notdir $(matrix)))_p$(p)
endif

ifeq ($(IMPL),pthread)
PSPGEMM_SRC = pthread.c
PSPGEMM_LIBS = -lpthread
else ifeq ($(IMPL),omp)
PSPGEMM_SRC = omp.c
PSPGEMM_LIBS = -fopenmp
else ifeq ($(IMPL),balance)
PSPGEMM_SRC = balance.c
PSPGEMM_LIBS = -lpthread
else ifeq ($(IMPL),balance_v2)
PSPGEMM_SRC = balance_v2.c
PSPGEMM_LIBS = -lpthread
else
PSPGEMM_SRC =
PSPGEMM_LIBS = 
endif

full: clean
	sbatch --partition=$(PARTITION) \
		--job-name=$(RUN_TAG) \
		--output=slurn_%j_$(RUN_TAG).out \
		--error=slurn_%j_$(RUN_TAG).err \
		--cpus-per-task=$(p) \
		--export=MODE="$(MODE)",IMPL="$(IMPL)",matrix="$(matrix)",p="$(RUN_P)",test="$(test)" run.sh

pspgemm: $(PSPGEMM_SRC)
	$(CC) $(CFLAGS) $(TEST_FLAG) -o $@ $(PSPGEMM_SRC) $(LDFLAGS) $(PSPGEMM_LIBS)

serial: serial.c
	$(CC) $(CFLAGS) -DSERIAL_STANDALONE -o $@ serial.c $(LDFLAGS)

clean:
	rm -f pspgemm serial

realclean: clean
	rm -f *.err *.out