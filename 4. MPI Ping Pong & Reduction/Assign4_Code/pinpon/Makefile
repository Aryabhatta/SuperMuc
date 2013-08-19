# Intel compiler
CC =  icc

CFLAGS = -g

MPICC = mpicc

all: PinPon

PinPon : pingpong.o
	$(MPICC) $(CFLAGS) -o $@ $+ 

%.o : %.c
	$(MPICC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o PinPon

remake : clean all
