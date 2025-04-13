#include "include/datapoint.h"

#include <vector>
#include <cmath>
#include <cstring> // for memcpy

// Constructor
DataPoint::DataPoint(double *values, int dim) {
    this->data = values; this->dim = dim;
}

// Destructor
DataPoint::~DataPoint() {
    delete[] data; // Free the allocated memory
}

// Accessor
double DataPoint::operator[](int index) const {
    return data[index];
}

// Copy method
void DataPoint::copy(double *arr, int n) const {
    std::memcpy(arr, data, n * sizeof(double));
}