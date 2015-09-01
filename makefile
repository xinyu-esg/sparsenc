# makefile for this mini-project
#

CC = gcc
#CFLAGS = -fdiagnostics-color=auto -std=c99 -g -lm
CFLAGS = -std=c99 -g -lm

TEST_GGDEC = test.GGdecoder.o bipartite.o gncEncoder.o gncGGDecoder.o galois.o gaussian.o common.o
TEST_OADEC = test.OAdecoder.o bipartite.o gncEncoder.o gncOADecoder.o galois.o gaussian.o common.o pivoting.o


test.GGdecoder: $(TEST_GGDEC)
	$(CC) -o $@ $^ $(CFLAGS)

test.OAdecoder: $(TEST_OADEC)
	$(CC) -o $@ $^ $(CFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)


.PHONY: clean

clean:
	rm -f *.o test.GGdecoder test.OAdecoder
