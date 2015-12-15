CUDA := g++
C_FLAGS := -lm -lfftw3

all: pic_leap

pic_leap: pic_leap.cpp
	$(CUDA) pic_leap.cpp -o pic_leap.mio $(C_FLAGS)

clean:
	rm -rf *.mio *.data *.out
