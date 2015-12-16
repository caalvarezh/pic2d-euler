CC := g++
C_FLAGS := -lm -lfftw3 -lOpenCL

all: pic_leap

pic_leap: pic_leap.cpp
	$(CC) pic_leap.cpp -o pic_leap.mio $(C_FLAGS)

clean:
	rm -rf *.mio *.data *.out
