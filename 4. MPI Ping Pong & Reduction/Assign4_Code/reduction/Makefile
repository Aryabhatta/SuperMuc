# Intel compiler
CC =  icc

CFLAGS = -g

MPICC = mpicc

all: Reduction

Reduction : reduction.o
	$(MPICC) $(CFLAGS) -o $@ $+ -lm

%.o : %.c
	$(MPICC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o Reduction

remake : clean all
