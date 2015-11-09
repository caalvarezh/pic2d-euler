CXX := g++
CUDA := nvcc
C_FLAGS := -lm -lfftw3
CUDA_FLAGS := -arch=compute_35 -lcufft

all: pic_leap pic_generator pic_test

parallel: pic_leap

sec: secuential

test: pic_test

pic_leap: pic_leap.cu
	$(CUDA) pic_leap.cu -o pic_leap.mio $(C_FLAGS) $(CUDA_FLAGS) 

secuential: pic_leap.cpp
	$(CXX) pic_leap.cpp -o pic_leap.mio $(C_FLAGS)

generator: pic_generator.cu
	$(CUDA) pic_generator.cu -o pic_generator.mio $(C_FLAGS) $(CUDA_FLAGS) 

pic_test: pic_test.cu
	$(CUDA) $(CUDA_FLAGS) pic_test.cu -o pic_test.mio $(C_FLAGS)

clean:
	rm -rf *.mio *.data *.out
