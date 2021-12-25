CC=gcc
CFLAGS=-Wall -Wextra
OFLAGS=-O3 -march=native -mtune=native
DFLAGS=-g
LFLAGS=

TARGET=main

.PHONY: all clean

all: $(TARGET)

$(TARGET): main.o common.o lennard_jones.o
	$(CC) $(LFLAGS) $^ -o $@

%.o: %.c
	$(CC) -c $(CFLAGS) $(OFLAGS) $(DFLAGS) $< -o $@

main.c: lennard_jones.c lennard_jones.h common.c common.h helper.h

lennard_jones.c: lennard_jones.h common.c common.h helper.h

common.c: common.h helper.h

clean:
	rm -Rf *~ *.o $(TARGET)
