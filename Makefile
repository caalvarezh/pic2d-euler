CXX := g++
CUDA := nvcpp
C_FLAGS := -lm -lfftw3
CUDA_FLAGS :=

all: pic_leap pic_generator pic_test

test: pic_test pic_generator

pic_leap: pic_leap.cpp
	$(CXX) $(C_FLAGS) pic_leap.cpp -o pic_leap.mio

pic_generator: pic_generator.cpp
	$(CXX) $(C_FLAGS) pic_generator.cpp -o pic_generator.mio

pic_test: pic_test.cpp
	$(CUDA) $(CUDA_FLAGS) $(C_FLAGS) pic_test.cu -o pic_test.mio

clean:
	rm -rf *.mio *.data
