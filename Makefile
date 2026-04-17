CC      = gcc
CFLAGS  = -O3 -Wall -Wextra -march=native
LDFLAGS = -lm
PARTITION = ece6750

MODE ?= single
IMPL ?= pthreads
matrix ?= small.mtx
p ?= 1
test ?= 0

ifeq ($(test),1)
TESTFLAG = -DTEST
else
TESTFLAG =
endif

.PHONY: all clean full single parallel

all:
	@echo "Use: make full MODE=single matrix=small.mtx"
	@echo " or: make full MODE=parallel IMPL=pthreads matrix=even.mtx p=8"

single: serial

parallel:
ifeq ($(IMPL),pthreads)
	$(MAKE) pthreads
else ifeq ($(IMPL),omp)
	$(MAKE) omp
else ifeq ($(IMPL),balance)
	$(MAKE) balance
else
	@echo "Unknown IMPL=$(IMPL). Use pthreads, omp, or balance."
	@exit 1
endif

serial: serial.c
	$(CC) $(CFLAGS) $(TESTFLAG) serial.c -o serial $(LDFLAGS)

pthreads: pthread.c
	$(CC) $(CFLAGS) $(TESTFLAG) pthread.c -o pthread $(LDFLAGS) -lpthread

omp: omp.c
	$(CC) $(CFLAGS) $(TESTFLAG) omp.c -o omp $(LDFLAGS) -fopenmp

balance: balance.c
	$(CC) $(CFLAGS) $(TESTFLAG) balance.c -o balance $(LDFLAGS) -lpthread

full:
	sbatch --partition=$(PARTITION) \
	       --cpus-per-task=$(p) \
	       --export=MODE=$(MODE),IMPL=$(IMPL),matrix=$(matrix),p=$(p),test=$(test) \
	       run.sh

clean:
	rm -f serial pthread omp balance
	rm -f slurm-*.out slurm-*.err