# Intel compiler
#CC =  icc

CFLAGS = -g -O2 -fargument-noalias 

MPICC = mpicc

all: heat

heat : heat.o input.o misc.o timing.o relax_gauss.o relax_jacobi.o
	$(MPICC) $(CFLAGS) -o $@ $+ -lm

#-L/lrz/sys/tools/papi/4.9.0-rapl/lib/libpapi.a -L/lrz/sys/tools/papi/4.9.0-rapl/lib/libpfm.a

%.o : %.c %.h
	$(MPICC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o heat *~ *.ppm

remake : clean all
