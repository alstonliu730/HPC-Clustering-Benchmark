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

# Targets
TARGETS = main dbscan-cpu dbscan-cuda

# main
#all:
#	make $(TARGETS)

main: main.cpp $(UTILS_DIR)/*.cpp $(SRC_DIR)/*.cpp
	$(CXX) $(CFLAGS) $(DEBUG) $(OMP_FLAGS) -I $(INC_DIR)/ -I $(UTILS_DIR)/ -o $@.exe $^

dbscan-cpu:$(SRC_DIR)/dbscan-cpu.cpp
	$(CXX) $(CFLAGS) $(OMP_FLAGS) -I $(INC_DIR)/ -o $@ $<

# dbscan-cuda:

clean:
	rm -f $(TARGETS)
