CC=gcc

CFLAGS=-Wall -g -O3 -msse4
LDFLAGS=-g -Wall -O3 -msse4
LIBFLAGS=-lm -ljpeg -lpng
LIBAFLAGS=-static /usr/lib64/libjpeg.a /usr/lib64/libpng.a /usr/lib64/libz.a /usr/lib64/libm.a

all: fastdeepflow

fastdeepflow: deepflow.o image.o io.o opticalflow_aux.o opticalflow.o solver.o
	$(CC) $(LDFLAGS) $(LIBFLAGS) -o $@ $^

fastdeepflow-static: deepflow.o image.o io.o opticalflow_aux.o opticalflow.o solver.o
	$(CC) $(LIBFLAGS) -o $@ $^ $(LIBAFLAGS)

%.o: %.c
	$(CC) -o $@ $(CFLAGS) -c $+ 

clean:
	rm -f *.o deepflow
