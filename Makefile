# Tested that both gcc and clang-12 will successfully compile all executables
#
# clang-12 produces code that is slightly faster than gcc, but must use the
# -fno-inline flag when compiling s25
#
# Use gcc for consistent optimization behavior

#all: a25 s25 v25 525
all: s25 v25

CC=clang-12
#CC=gcc

CFLAGS=-O3 -march=native -Wall
LIBS=-lpthread

#a25: a25.c utilities.h Makefile
#	$(CC) -fno-inline $(CFLAGS) -o $@ a25.c $(LIBS)

s25: s25.c utilities.h Makefile
	$(CC) $(CFLAGS) -o $@ s25.c $(LIBS)

v25: v25.c utilities.h Makefile
	$(CC) $(CFLAGS) -o $@ v25.c $(LIBS)

#525: 525.c utilities.h Makefile
#	$(CC) $(CFLAGS) -o $@ 525.c $(LIBS)

check:
	/bin/sh ./check.sh
