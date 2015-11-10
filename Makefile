######################################################
# Makefile for libgnc
# Ye Li
# leeyee.seu@gmail.com
######################################################

TOP = .
SRCDIR := src
OBJDIR := src
INCLUDEDIR = include
INC_PARMS = $(INCLUDEDIR:%=-I%)

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
	SED = gsed
	CC  = gcc-5
endif
ifeq ($(UNAME), Linux)
	SED = sed
	CC  = gcc
endif

#CC = gcc
CFLAGS0 = -Winline -std=c99 -lm
CFLAGS1 = -O3 -DNDEBUG $(INC_PARMS)  -mssse3 -DINTEL_SSSE3
#CFLAGS2 = -lm
#CFLAGS = -O3 -I$(INCLUDEDIR) -mssse3 -DINTEL_SSSE3
#CFLAGS = -std=c99 -g -lm

vpath %.h src include
vpath %.c src examples

DEFS   := common.h bipartite.h gncEncoder.h galois.h
GNCENC  := $(OBJDIR)/common.o $(OBJDIR)/bipartite.o $(OBJDIR)/gncEncoder.o $(OBJDIR)/galois.o $(OBJDIR)/gaussian.o
RECODER := $(OBJDIR)/gncRecoder.o 
GGDEC   := $(OBJDIR)/gncGGDecoder.o 
OADEC   := $(OBJDIR)/gncOADecoder.o $(OBJDIR)/pivoting.o
BDDEC   := $(OBJDIR)/gncBandDecoder.o $(OBJDIR)/pivoting.o
CBDDEC  := $(OBJDIR)/gncCBDDecoder.o
DECDEFS := gncGGDecoder.h gncOADecoder.h gncBandDecoder.h gncCBDDecoder.h

.PHONY: all
all: band.OA.example band.GG.example band.BD.example band.CBD.example rand.GG.example rand.OA.example

libgnc.so: $(GNCENC) $(GGDEC) $(OADEC) $(BDDEC) $(CBDDEC) $(RECODER)
	$(CC) -shared -o libgnc.so $^
	
#GGband.example
band.GG.example: libgnc.so test.GGdecoder.c
	$(SED) -i 's/[^ ]*GNC_CODE/BAND_GNC_CODE/' examples/test.GGdecoder.c
	$(CC) -L. -lgnc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#GGrand.example
rand.GG.example: libgnc.so test.GGdecoder.c
	$(SED) -i 's/[^ ]*GNC_CODE/RAND_GNC_CODE/' examples/test.GGdecoder.c
	$(CC) -L. -lgnc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#OAband.example
band.OA.example: libgnc.so test.OAdecoder.c 
	$(SED) -i 's/[^ ]*GNC_CODE/BAND_GNC_CODE/' examples/test.OAdecoder.c
	$(CC) -L. -lgnc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#OAband statically linked, not using libgnc.so
band.OA.static: $(GNCENC) $(OADEC) examples/test.OAdecoder.c
	$(SED) -i 's/[^ ]*GNC_CODE/BAND_GNC_CODE/' examples/test.OAdecoder.c
	$(CC) -o $@ $(CFLAGS0) $(CFLAGS1) $^
#OArand.example
rand.OA.example: libgnc.so test.OAdecoder.c 
	$(SED) -i 's/[^ ]*GNC_CODE/RAND_GNC_CODE/' examples/test.OAdecoder.c
	$(CC) -L. -lgnc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#OAwind.example
wind.OA.example: libgnc.so test.OAdecoder.c 
	$(SED) -i 's/[^ ]*GNC_CODE/WINDWRAP_GNC_CODE/' examples/test.OAdecoder.c
	$(CC) -L. -lgnc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#BDband.example
band.BD.example: libgnc.so test.bandDecoder.c 
	$(SED) -i 's/[^ ]*GNC_CODE/BAND_GNC_CODE/' examples/test.bandDecoder.c
	$(CC) -L. -lgnc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#BDband statically linked, not using libgnc.so
band.BD.static: $(GNCENC) $(BDDEC) examples/test.bandDecoder.c
	$(SED) -i 's/[^ ]*GNC_CODE/BAND_GNC_CODE/' examples/test.bandDecoder.c
	$(CC) -o $@ $(CFLAGS0) $(CFLAGS1) $^
#CBDband.example
band.CBD.example: libgnc.so test.CBDDecoder.c 
	$(SED) -i 's/[^ ]*GNC_CODE/BAND_GNC_CODE/' examples/test.CBDDecoder.c
	$(CC) -L. -lgnc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#Recoder with band code and CBD decoder, TRIV_SCHED
recoder.CBD.trivSched: libgnc.so test.2hopRecoder.CBD.c
	$(SED) -i 's/[^ ]*GNC_CODE/BAND_GNC_CODE/' examples/test.2hopRecoder.CBD.c
	$(SED) -i 's/[^ ]*_SCHED/TRIV_SCHED/' examples/test.2hopRecoder.CBD.c
	$(CC) -L. -lgnc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#Recoder with band code and CBD decoder, RAND_SCHED
recoder.CBD.randSched: libgnc.so test.2hopRecoder.CBD.c
	$(SED) -i 's/[^ ]*GNC_CODE/BAND_GNC_CODE/' examples/test.2hopRecoder.CBD.c
	$(SED) -i 's/[^ ]*_SCHED/RAND_SCHED/' examples/test.2hopRecoder.CBD.c
	$(CC) -L. -lgnc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#Recoder with band code and CBD decoder, MLPI_SCHED
recoder.CBD.mlpiSched: libgnc.so test.2hopRecoder.CBD.c
	$(SED) -i 's/[^ ]*GNC_CODE/BAND_GNC_CODE/' examples/test.2hopRecoder.CBD.c
	$(SED) -i 's/[^ ]*_SCHED/MLPI_SCHED/' examples/test.2hopRecoder.CBD.c
	$(CC) -L. -lgnc -o $@ $(CFLAGS0) $(CFLAGS1) $^

$(OBJDIR)/%.o: $(OBJDIR)/%.c $(DEFS) $(GGDEFS)
	$(CC) -c -fpic -o $@ $< $(CFLAGS0) $(CFLAGS1)

.PHONY: clean
clean:
	rm -f *.o $(OBJDIR)/*.o *.example *.static recoder.CBD.randSched recoder.CBD.trivSched recoder.CBD.mlpiSched libgnc.so

.PHONY: install
install:
	cp libgnc.so /usr/lib/

.PHONY: uninstall
uninstall:
	rm -f /usr/lib/libgnc.so
