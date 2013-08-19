# Intel compiler
CC =  icc

CFLAGS = -g -O2 -openmp $$PAPI_INC

MPICC = mpicc

all: heat

heat : heat.o input.o misc.o timing.o relax_gauss.o relax_jacobi.o
	$(CC) $(CFLAGS) -o $@ $+ -lm $$PAPI_LIB
#-L/lrz/sys/tools/papi/4.9.0-rapl/lib/libpapi.a -L/lrz/sys/tools/papi/4.9.0-rapl/lib/libpfm.a



%.o : %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o heat *~ *.ppm

remake : clean all
