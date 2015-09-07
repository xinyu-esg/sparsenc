# makefile for this mini-project
#

CC = gcc
#CFLAGS = -fdiagnostics-color=auto -std=c99 -g -lm
CFLAGS = -std=c99 -g -lm

TEST_GGDEC = test.GGdecoder.o bipartite.o gncEncoder.o gncGGDecoder.o galois.o gaussian.o common.o
TEST_GGDEC_FILE = test.GGdecoder.file.o bipartite.o gncEncoder.o gncGGDecoder.o galois.o gaussian.o common.o
TEST_OADEC = test.OAdecoder.o bipartite.o gncEncoder.o gncOADecoder.o galois.o gaussian.o common.o pivoting.o
TEST_OADEC_FILE = test.OAdecoder.file.o bipartite.o gncEncoder.o gncOADecoder.o galois.o gaussian.o common.o pivoting.o
TEST_BDDEC = test.bandDecoder.o bipartite.o gncEncoder.o gncBandDecoder.o galois.o gaussian.o common.o pivoting.o
TEST_BDDEC_FILE = test.bandDecoder.file.o bipartite.o gncEncoder.o gncBandDecoder.o galois.o gaussian.o common.o pivoting.o


test.GGdecoder: $(TEST_GGDEC)
	$(CC) -o $@ $^ $(CFLAGS)

test.GGdecoder.file: $(TEST_GGDEC_FILE)
	$(CC) -o $@ $^ $(CFLAGS)

test.OAdecoder: $(TEST_OADEC)
	$(CC) -o $@ $^ $(CFLAGS)

test.OAdecoder.file: $(TEST_OADEC_FILE)
	$(CC) -o $@ $^ $(CFLAGS)

test.bandDecoder: $(TEST_BDDEC)
	$(CC) -o $@ $^ $(CFLAGS)

test.bandDecoder.file: $(TEST_BDDEC_FILE)
	$(CC) -o $@ $^ $(CFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)


.PHONY: clean

clean:
	rm -f *.o test.GGdecoder test.GGdecoder.file test.OAdecoder test.OAdecoder.file test.bandDecoder test.bandDecoder.file
