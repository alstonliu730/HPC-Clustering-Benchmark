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
UTILS_DIR = utils/
TARGETS = main dbscan-cpu dbscan-cuda

# main
all: 
	make $(TARGETS)
main: $(SRC_DIR)main.cpp
    $(CXX) $(CFLAGS) $(OMP_FLAGS) -I$(INC_DIR) -o $@ $<
	
dbscan-cpu:$(SRC_DIR)dbscan-cpu.cpp
	$(CXX) $(CFLAGS) $(OMP_FLAGS) -I $(INC_DIR)*.h -o $@ $<

# dbscan-cuda:

clean:
	rm -f $(TARGETS)
