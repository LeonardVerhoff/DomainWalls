CC = cc #gcc

CFLAGS = -g #-O3		# debugging information, optimization level
CPPFLAGS = -I.		# Flags to the C preprocessor, here include files from the current directory
LOADLIBES = -lm		# Libraries to load, here libm - the math library
VPATH = ../src		# list of directories that make should search

all:
	$(CC) -g -Wall monte_carlo.c ran.c common.c -o monte_carlo -lm 



