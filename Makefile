CC = gcc
CFLAGS = -Wall -Wextra -O2
CXX = g++

# CUDA compiler
NVCC = nvcc
CUDA_FLAGS = -O3
OMP_FLAGS = -fopenmp

# Source directory
SRC_DIR = src/
INC_DIR = include/
TARGETS = main dbscan-cpu dbscan-cuda

# main
all: 
	make $(TARGETS)

dbscan-cpu: $(INC_DIR)dbscan.h $(SRC_DIR)dbscan-cpu.cpp 
	$(CXX) $(CFLAGS) -I$(INC_DIR) -o $@ $<

dbscan-cuda:

clean:
	rm -f $(TARGETS)
