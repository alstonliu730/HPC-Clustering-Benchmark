CC = gcc
CFLAGS = -Wall -Wextra -O2
CXX = g++
DEBUG = -g

# CUDA compiler
NVCC = nvcc
CUDA_FLAGS = -std=c++11 -O3
OMP_FLAGS = -fopenmp

# Source directory
SRC_DIR = src
INC_DIR = include
UTILS_DIR = utils
FLANN_DIR = flann/src/cpp
LZ4_DIR = external/lz4/lib
LZ4_LIB = $(LZ4_DIR)/liblz4.a

# Targets
TARGETS = main dbscan-cpu dbscan-cuda

# main
#all:
#	make $(TARGETS)

main: main.cpp $(UTILS_DIR)/*.cpp $(SRC_DIR)/*.cpp
	$(CXX) $(CFLAGS) $(DEBUG) $(OMP_FLAGS) \
	-I $(INC_DIR)/ -I $(UTILS_DIR)/ -I $(FLANN_DIR)/ -I $(LZ4_DIR) \
	-o $@.exe $^ $(LZ4_LIB)


dbscan-cuda: dbscan.cu
	$(NVCC) $(CUDA_FLAGS) -o $@ 

clean:
	rm -f $(TARGETS)
