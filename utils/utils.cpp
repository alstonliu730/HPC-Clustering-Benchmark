#include "utils.h"
#include <iostream>
#include <fstream>
#include <cstddef> // for size_t

size_t import_data(char* filename, std::vector<DataPoint>& points) {
    FILE* file = fopen(filename, "r");

    if (file == NULL) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return 0; // Return 0 if file cannot be opened
    }

    // read file until EOF
    double x, y, z;
    size_t count = 0;

    while (fscanf(file, "%lf, %lf, %lf", &x, &y, &z) == 3) {
        double* values = new double[3]; // Allocate memory for 3D point
        values[0] = x;
        values[1] = y;
        values[2] = z;
        points.push_back(DataPoint(values, 3)); // Add point to vector
        count++;
    }
}