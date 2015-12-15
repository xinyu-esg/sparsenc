######################################################
# Makefile for libsnc
# Ye Li
# leeyee.seu@gmail.com
######################################################

TOP = .
SRCDIR := src
OBJDIR := src
INCLUDEDIR = include src
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

DEFS    := snc.h common.h galois.h decoderGG.h decoderOA.h decoderBD.h decoderCBD.h
GNCENC  := $(OBJDIR)/common.o $(OBJDIR)/bipartite.o $(OBJDIR)/sncEncoder.o $(OBJDIR)/galois.o $(OBJDIR)/gaussian.o
RECODER := $(OBJDIR)/sncRecoder.o 
DECODER := $(OBJDIR)/sncDecoder.o
GGDEC   := $(OBJDIR)/decoderGG.o 
OADEC   := $(OBJDIR)/decoderOA.o $(OBJDIR)/pivoting.o
BDDEC   := $(OBJDIR)/decoderBD.o $(OBJDIR)/pivoting.o
CBDDEC  := $(OBJDIR)/decoderCBD.o

.PHONY: all
all: sncDecoder sncDecoderFile sncRecoder2Hop

libsnc.so: $(GNCENC) $(GGDEC) $(OADEC) $(BDDEC) $(CBDDEC) $(RECODER) $(DECODER)
	$(CC) -shared -o libsnc.so $^
	
#Test snc decoder
sncDecoder: libsnc.so test.decoders.c 
	$(CC) -L. -lsnc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#Test decoder for files
sncDecoderFile: libsnc.so test.file.decoders.c 
	$(CC) -L. -lsnc -o $@ $(CFLAGS0) $(CFLAGS1) $^
#Test recoder
sncRecoder2Hop: libsnc.so test.2hopRecoder.c
	$(CC) -L. -lsnc -o $@ $(CFLAGS0) $(CFLAGS1) $^

$(OBJDIR)/%.o: $(OBJDIR)/%.c $(DEFS)
	$(CC) -c -fpic -o $@ $< $(CFLAGS0) $(CFLAGS1)

.PHONY: clean
clean:
	rm -f *.o $(OBJDIR)/*.o libsnc.so sncDecoder sncDecoderFile sncRecoder2Hop

install: libsnc.so
	cp include/snc.h /usr/include/
	if [[ `uname -i` == "x86_64" ]]; then \
		cp libsnc.so /usr/lib64/; \
	else \
		cp libsnc.so /usr/lib/; \
	fi

.PHONY: uninstall
uninstall:
	rm -f /usr/include/snc.h
	if [[ `uname -i` == "X86_64" ]]; then \
		rm -f /usr/lib64/libsnc.so; \
	else \
		rm -f /usr/lib/libsnc.so; \
	fi
