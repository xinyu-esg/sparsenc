######################################################
# Makefile for libslnc
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

DEFS   := common.h bipartite.h slncEncoder.h galois.h
GNCENC  := $(OBJDIR)/common.o $(OBJDIR)/bipartite.o $(OBJDIR)/slncEncoder.o $(OBJDIR)/galois.o $(OBJDIR)/gaussian.o
RECODER := $(OBJDIR)/slncRecoder.o 
GGDEC   := $(OBJDIR)/slncGGDecoder.o 
OADEC   := $(OBJDIR)/slncOADecoder.o $(OBJDIR)/pivoting.o
BDDEC   := $(OBJDIR)/slncBandDecoder.o $(OBJDIR)/pivoting.o
CBDDEC  := $(OBJDIR)/slncCBDDecoder.o
DECDEFS := slncGGDecoder.h slncOADecoder.h slncBandDecoder.h slncCBDDecoder.h

.PHONY: all
all: band.OA.example band.GG.example band.BD.example band.CBD.example rand.GG.example rand.OA.example

libslnc.so: $(GNCENC) $(GGDEC) $(OADEC) $(BDDEC) $(CBDDEC) $(RECODER)
	$(CC) -shared -o libslnc.so $^
	
#GGband.example
band.GG.example: libslnc.so test.GGdecoder.c
	$(SED) -i 's/[^ ]*GNC_CODE/BAND_GNC_CODE/' examples/test.GGdecoder.c
	$(CC) -L. -lslnc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#GGrand.example
rand.GG.example: libslnc.so test.GGdecoder.c
	$(SED) -i 's/[^ ]*GNC_CODE/RAND_GNC_CODE/' examples/test.GGdecoder.c
	$(CC) -L. -lslnc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#OAband.example
band.OA.example: libslnc.so test.OAdecoder.c 
	$(SED) -i 's/[^ ]*GNC_CODE/BAND_GNC_CODE/' examples/test.OAdecoder.c
	$(CC) -L. -lslnc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#OAband statically linked, not using libslnc.so
band.OA.static: $(GNCENC) $(OADEC) examples/test.OAdecoder.c
	$(SED) -i 's/[^ ]*GNC_CODE/BAND_GNC_CODE/' examples/test.OAdecoder.c
	$(CC) -o $@ $(CFLAGS0) $(CFLAGS1) $^
#OArand.example
rand.OA.example: libslnc.so test.OAdecoder.c 
	$(SED) -i 's/[^ ]*GNC_CODE/RAND_GNC_CODE/' examples/test.OAdecoder.c
	$(CC) -L. -lslnc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#OAwind.example
wind.OA.example: libslnc.so test.OAdecoder.c 
	$(SED) -i 's/[^ ]*GNC_CODE/WINDWRAP_GNC_CODE/' examples/test.OAdecoder.c
	$(CC) -L. -lslnc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#BDband.example
band.BD.example: libslnc.so test.bandDecoder.c 
	$(SED) -i 's/[^ ]*GNC_CODE/BAND_GNC_CODE/' examples/test.bandDecoder.c
	$(CC) -L. -lslnc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#BDband statically linked, not using libslnc.so
band.BD.static: $(GNCENC) $(BDDEC) examples/test.bandDecoder.c
	$(SED) -i 's/[^ ]*GNC_CODE/BAND_GNC_CODE/' examples/test.bandDecoder.c
	$(CC) -o $@ $(CFLAGS0) $(CFLAGS1) $^
#CBDband.example
band.CBD.example: libslnc.so test.CBDDecoder.c 
	$(SED) -i 's/[^ ]*GNC_CODE/BAND_GNC_CODE/' examples/test.CBDDecoder.c
	$(CC) -L. -lslnc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#Recoder with band code and CBD decoder, TRIV_SCHED
recoder.CBD.trivSched: libslnc.so test.2hopRecoder.CBD.c
	$(SED) -i 's/[^ ]*GNC_CODE/BAND_GNC_CODE/' examples/test.2hopRecoder.CBD.c
	$(SED) -i 's/[^ ]*_SCHED/TRIV_SCHED/' examples/test.2hopRecoder.CBD.c
	$(CC) -L. -lslnc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#Recoder with band code and CBD decoder, RAND_SCHED
recoder.CBD.randSched: libslnc.so test.2hopRecoder.CBD.c
	$(SED) -i 's/[^ ]*GNC_CODE/BAND_GNC_CODE/' examples/test.2hopRecoder.CBD.c
	$(SED) -i 's/[^ ]*_SCHED/RAND_SCHED/' examples/test.2hopRecoder.CBD.c
	$(CC) -L. -lslnc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#Recoder with band code and CBD decoder, MLPI_SCHED
recoder.CBD.mlpiSched: libslnc.so test.2hopRecoder.CBD.c
	$(SED) -i 's/[^ ]*GNC_CODE/BAND_GNC_CODE/' examples/test.2hopRecoder.CBD.c
	$(SED) -i 's/[^ ]*_SCHED/MLPI_SCHED/' examples/test.2hopRecoder.CBD.c
	$(CC) -L. -lslnc -o $@ $(CFLAGS0) $(CFLAGS1) $^

$(OBJDIR)/%.o: $(OBJDIR)/%.c $(DEFS) $(GGDEFS)
	$(CC) -c -fpic -o $@ $< $(CFLAGS0) $(CFLAGS1)

.PHONY: clean
clean:
	rm -f *.o $(OBJDIR)/*.o *.example *.static recoder.CBD.randSched recoder.CBD.trivSched recoder.CBD.mlpiSched libslnc.so

.PHONY: install
install:
	cp libslnc.so /usr/lib/

.PHONY: uninstall
uninstall:
	rm -f /usr/lib/libslnc.so
