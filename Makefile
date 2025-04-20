# Compiler and flags
CXX = g++
MPI_CXX = mpicxx
CFLAGS = -Wall -Wextra -O2 -std=c++11
LDFLAGS = -fopenmp

# MPI flag (set to 1 to enable MPI, 0 to disable)
USE_MPI = 0

# Source directories
SRC_DIR = src
INC_DIR = include
UTILS_DIR = utils
FLANN_DIR = flann/src/cpp
LZ4_DIR = external/lz4/lib
LZ4_LIB = $(LZ4_DIR)/liblz4.a

# Object files
OBJS = $(SRC_DIR)/dbscan_cpu.cpp $(SRC_DIR)/datapoint.cpp $(SRC_DIR)/main.cpp
MPI_OBJS = $(SRC_DIR)/dbscan_cpu.cpp $(SRC_DIR)/datapoint.cpp main_mpi.cpp

# Include directories
INCLUDES = -I$(INC_DIR)/ -I$(UTILS_DIR)/ -I$(FLANN_DIR)/ -I$(LZ4_DIR)/

# MPI-specific flags and libraries
MPI_CFLAGS = -I/usr/include/mpi
MPI_LDFLAGS = -lmpi

# Target executable
TARGET = main 

# Build rules
ifeq ($(USE_MPI), 1)
	CXX = $(MPI_CXX)
	OBJS = $(MPI_OBJS)
    CFLAGS += $(MPI_CFLAGS)
    LDFLAGS += $(MPI_LDFLAGS)
    DEFINES += -DUSE_MPI
endif

all: $(TARGET)

$(TARGET): $(OBJS)
    $(CXX) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS)

clean:
    rm -f $(TARGET) *.o
