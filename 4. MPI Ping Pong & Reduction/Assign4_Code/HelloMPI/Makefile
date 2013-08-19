# Intel compiler
CC =  icc

CFLAGS = -g

MPICC = mpicc

all: Hello

Hello : HelloMpi.o
	$(MPICC) $(CFLAGS) -o $@ $+ 

%.o : %.c
	$(MPICC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o Hello

remake : clean all
