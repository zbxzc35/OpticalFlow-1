CC=gcc

CFLAGS=-Wall -g -O3 -msse4
LDFLAGS=-g -Wall -O3 -msse4
LIBFLAGS=-lm -ljpeg -lpng
#LIBAFLAGS=-static /usr/lib/x86_64-linux-gnu/libjpeg.a /usr/lib/x86_64-linux-gnu/libpng.a /usr/lib/x86_64-linux-gnu/libz.a /usr/lib/x86_64-linux-gnu/libm.a

all:fastdeepflow
#all: fastdeepflow-static

fastdeepflow: deepflow.o image.o io.o opticalflow_aux.o opticalflow.o solver.o
	$(CC) $(LDFLAGS) $(LIBFLAGS) -o $@ $^ $(LIBFLAGS)

#fastdeepflow-static: deepflow.o image.o io.o opticalflow_aux.o opticalflow.o solver.o
#	$(CC) $(LIBFLAGS) -o $@ $^ $(LIBAFLAGS)

%.o: %.c
	$(CC) -o $@ $(CFLAGS) -c $+ 

clean:
	rm -f *.o deepflow
