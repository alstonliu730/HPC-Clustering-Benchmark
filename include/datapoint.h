#ifndef _DATAPOINT_H_
#define _DATAPOINT_H_

#include <cstring> // for memcpy

class DataPoint {
    private:
        float *data; // Pointer to the data array
        int dim; // Dimension of the data point
        bool visited = false; // Flag to indicate if the point has been visited
    public:
        // Constructor
        DataPoint(float *values, int dim);

        // Copy constructor
        DataPoint(const DataPoint& other);
        
        // Destructor
        ~DataPoint();

        // Accessor
        float operator[](int index) const;

        // Copy assignment operator
        DataPoint& operator=(const DataPoint& other);

        // Copy method
        void copy(float *arr, int n) const;

        // Get Dimension
        int get_dim() const { return dim; }
        
        // Get Data
        float* get_data() const { return data; }
};

#endif // DATAPOINT_H_