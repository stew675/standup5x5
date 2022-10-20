# Tested that both gcc and clang-12 will successfully compile all executables
#
# clang-12 produces code that is slightly faster than gcc, but must use the
# -fno-inline flag when compiling s25
#
# Use gcc for consistent optimization behavior

all: a25 s25 v25 525

CC?=clang
CXX?=clang++
#CC=gcc
#CXX=g++

CFLAGS?=-O3 -march=native -Wall
CXXFLAGS?=-std=c++11 -O3 -Wall
LIBS=-lpthread

build_header: build_header.cpp Makefile
	$(CXX) $(CXXFLAGS) -o $@ build_header.cpp

words_alpha.h: build_header words_alpha.txt Makefile
	./build_header words_alpha.txt

words_alpha_five.h: build_header words_alpha_five.txt Makefile
	./build_header words_alpha_five.txt

nyt_wordle.h: build_header nyt_wordle.txt Makefile
	./build_header nyt_wordle.txt

a25: a25.c utilities.h words_alpha.h words_alpha_five.h nyt_wordle.h Makefile
	$(CC) $(CFLAGS) -o $@ a25.c $(LIBS)

s25: s25.c utilities.h words_alpha.h words_alpha_five.h nyt_wordle.h Makefile
	$(CC) $(CFLAGS) -o $@ s25.c $(LIBS)

v25: v25.c utilities.h words_alpha.h words_alpha_five.h nyt_wordle.h Makefile
	$(CC) $(CFLAGS) -o $@ v25.c $(LIBS)

525: 525.c utilities.h words_alpha.h words_alpha_five.h nyt_wordle.h Makefile
	$(CC) $(CFLAGS) -o $@ 525.c $(LIBS)

check:
	/bin/sh ./check.sh
