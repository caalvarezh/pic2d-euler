CXX := g++
CUDA := nvcc
C_FLAGS := -lm -lfftw3
CUDA_FLAGS := -arch=compute_35

all: pic_leap pic_generator pic_test

test: pic_test

pic_leap: pic_leap.cpp
	$(CXX) pic_leap.cpp -o pic_leap.mio $(C_FLAGS)

pic_generator: pic_generator.cpp
	$(CXX)  pic_generator.cpp -o pic_generator.mio $(C_FLAGS)

pic_test: pic_test.cu
	$(CUDA) $(CUDA_FLAGS) pic_test.cu -o pic_test.mio $(C_FLAGS)

clean:
	rm -rf *.mio *.data *.out
