#ifndef _DATAPOINT_H_
#define _DATAPOINT_H_

#include <cstring> // for memcpy

class DataPoint {
    private:
        double *data; // Pointer to the data array
        int dim; // Dimension of the data point
    public:
        // Constructor
        DataPoint(double *values, int dim);
        
        // Destructor
        ~DataPoint();

        // Accessor
        double operator[](int index) const;

        // Copy method
        void copy(double *arr, int n) const;

        // Get Dimension
        int get_dim() const { return dim; }
        
        // Get Data
        double* get_data() const { return data; }
};

#endif // DATAPOINT_H_