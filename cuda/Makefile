CUDA := nvcc
C_FLAGS := -lm -lfftw3
CUDA_FLAGS := -arch=compute_35 -lcufft

all: pic_leap

pic_leap: pic_leap.cu
	$(CUDA) pic_leap.cu -o pic_leap.mio $(C_FLAGS) $(CUDA_FLAGS)

clean:
	rm -rf *.mio *.data *.out
