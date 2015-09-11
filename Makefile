######################################################
# Makefile for libgnc
# Ye Li
# leeyee.seu@gmail.com
######################################################

TOP = .
OBJDIR := src
INCLUDEDIR = include
SED = sed


CC = gcc
CFLAGS0 = -Winline -std=c99 -lm
CFLAGS = -std=c99 -I$(INCLUDEDIR) -lm 
#CFLAGS = -O3 -I$(INCLUDEDIR) -mssse3 -DINTEL_SSSE3
#CFLAGS = -std=c99 -g -lm

vpath %.h include
vpath %.c src examples

DEFS   := common.h bipartite.h gncEncoder.h galois.h
#GNCENC := $(addprefix, $(OBJDIR)/,common.o bipartite.o gncEncoder.o galois.o gaussian.o)
GNCENC  := $(OBJDIR)/common.o $(OBJDIR)/bipartite.o $(OBJDIR)/gncEncoder.o $(OBJDIR)/galois.o $(OBJDIR)/gaussian.o
GGDEC   := $(OBJDIR)/gncGGDecoder.o 
OADEC   := $(OBJDIR)/gncOADecoder.o $(OBJDIR)/pivoting.o
BDDEC   := $(OBJDIR)/gncBandDecoder.o $(OBJDIR)/pivoting.o
DECDEFS := gncGGDecoder.h gncOADecoder.h gncBandDecoder.h

.PHONY: all
all: band.OA.example band.GG.example band.BD.example rand.GG.example rand.OA.example
	
#GGband.example
band.GG.example: $(GNCENC) $(GGDEC) test.GGdecoder.c example_utils.c
	$(SED) -i 's/gnc_type\s=.*/gnc_type\ =\ BAND_GNC_CODE;/' examples/test.GGdecoder.c
	$(CC) -o $@ $(CFLAGS0) $(CFLAGS) $^
#GGrand.example
rand.GG.example: $(GNCENC) $(GGDEC) test.GGdecoder.c example_utils.c
	$(SED) -i 's/gnc_type\s=.*/gnc_type\ =\ RAND_GNC_CODE;/' examples/test.GGdecoder.c
	$(CC) -o $@ $(CFLAGS0) $(CFLAGS) $^
#OAband.example
band.OA.example: $(GNCENC) $(OADEC) test.OAdecoder.c example_utils.c
	$(SED) -i 's/gnc_type\s=.*/gnc_type\ =\ BAND_GNC_CODE;/' examples/test.OAdecoder.c
	$(CC) -o $@ $(CFLAGS0) $(CFLAGS) $^
#OArand.example
rand.OA.example: $(GNCENC) $(OADEC) test.OAdecoder.c example_utils.c
	$(SED) -i 's/gnc_type\s=.*/gnc_type\ =\ RAND_GNC_CODE;/' examples/test.OAdecoder.c
	$(CC) -o $@ $(CFLAGS0) $(CFLAGS) $^
#BDband.example
band.BD.example: $(GNCENC) $(BDDEC) test.bandDecoder.c example_utils.c
	$(SED) -i 's/gnc_type\s=.*/gnc_type\ =\ BAND_GNC_CODE;/' examples/test.bandDecoder.c
	$(CC) -o $@ $(CFLAGS0) $(CFLAGS) $^

$(OBJDIR)/%.o: $(OBJDIR)/%.c $(DEFS) $(GGDEFS)
	$(CC) -c -o $@ $< $(CFLAGS0) $(CFLAGS)

.PHONY: clean

clean:
	rm -f *.o $(OBJDIR)/*.o *.example
