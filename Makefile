#
# Makefile for Ens_Reweight_wExpdatda program
#

HOME=../
SRC=.
BIN=../bin
LIBDIR=-L/home/yamamori/work/programs/ABAMD_2014_05_01/lib #-L~/mylib -L~/lib
INCDIR=-I/home/yamamori/work/programs/ABAMD_2014_05_01/include #-I~/myinclude -I~/include

# definitions
#CC    = icc -pg -O3 -tpp6 -ipo
#CC    = icc -check=conversion,stack,uninit
#CC    = icc
CC    = gcc -g
#CC     = gcc

CFLAG = #-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -m64 -O2

OBJS =  main.o EF.o

LIBS =  EF.h

EXLIBS=-lm -lgc

TARGET =WHAM_MLE_oneD

.c.o:
	$(CC) $(INCDIR) -c $(CFLAG) $<;

all: WHAM_MLE_oneD install

# rules of generations
WHAM_MLE_oneD:  $(OBJS) $(LIBS)
	$(CC) $(CFLAG) -o $@ $(OBJS) $(LIBDIR) $(INCDIR) $(EXLIBS) ;

install: 
	cp $(TARGET) $(BIN) ; 

main.o:	EF.h
EF.o: EF.h 

clean: 
	rm $(OBJS); \
	rm $(TARGET);
