CC = gcc
CFLAGS = -Wall -Wextra -O2
CXX = g++
DEBUG = -g

# CUDA compiler
NVCC = nvcc
CUDA_FLAGS = -O3
OMP_FLAGS = -fopenmp

# Source directory
SRC_DIR = src
INC_DIR = include
UTILS_DIR = utils
FLANN_DIR = flann/src/cpp/flann

# Targets
TARGETS = main dbscan-cpu dbscan-cuda

# main
#all:
#	make $(TARGETS)

main: main.cpp $(UTILS_DIR)/*.cpp $(SRC_DIR)/*.cpp
	$(CXX) $(CFLAGS) $(DEBUG) $(OMP_FLAGS) -I $(INC_DIR)/\
	-I $(UTILS_DIR)/ -I $(FLANN_DIR)/ -o $@.exe $^

# dbscan-cuda:

clean:
	rm -f $(TARGETS)
