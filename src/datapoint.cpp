#include "datapoint.h"

#include <vector>
#include <cmath>
#include <cstring> // for memcpy

// Constructor
DataPoint::DataPoint(double *values, int dim) {
    this->data = new double[dim]; 
    this->dim = dim;
    std::memcpy(this->data, values, dim * sizeof(double)); // copy the values into the data array
}

// Copy constructor
DataPoint::DataPoint(const DataPoint& other) {
    this->dim = other.dim;
    this->data = new double[dim];
    std::memcpy(this->data, other.data, dim * sizeof(double));
}

// Copy assignment
DataPoint& DataPoint::operator=(const DataPoint& other) {
    if (this == &other) return *this; // self-assignment check

    delete[] this->data; // free existing memory
    this->dim = other.dim;
    this->data = new double[dim];
    std::memcpy(this->data, other.data, dim * sizeof(double));
    return *this;
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
