#
# Makefile for Ens_Reweight_wExpdatda program
#

HOME=. #../
SRC=.
BIN=. #../bin
LIBDIR=-L/Users/yamamoriyuu//onedrive/Research/work/programs/C/liblbfgs/lib
INCDIR=-I/Users/yamamoriyuu//onedrive/Research/work/programs/C/liblbfgs/include

# definitions
#CC    = icc -pg -O3 -tpp6 -ipo
#CC    = icc -check=conversion,stack,uninit
#CC    = icc
CC    = gcc -g
#CC     = gcc

CFLAG = #-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -m64 -O2

OBJS =  main.o error.o #func.o

LIBS =  error.h #func.h

EXLIBS = -llbfgs -lm -lgc

TARGET = ensRefEx

.c.o:
	$(CC) $(INCDIR) -c $(CFLAG) $<;

all: ensRefEx #install

# rules of generations

ensRefEx:  $(OBJS) $(LIBS)
	$(CC) $(CFLAG) -o $@ $(OBJS) $(LIBDIR) $(INCDIR) $(EXLIBS) ;

install: 
	cp $(TARGET) $(BIN) ; 

main.o: error.h #func.h
#func.o: func.h error.h
error.o: error.h 

clean: 
	rm $(OBJS); \
	rm $(TARGET);




