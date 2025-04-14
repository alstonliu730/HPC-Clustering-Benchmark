#include "utils.h"
#include "datapoint.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

size_t import_data(char* input_file, std::vector<DataPoint>& points) {
    std::ifstream file(input_file);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << input_file << std::endl;
        return 0;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Split the line by commas
        std::string delimiter = ",";
        size_t pos = 0;

        // create a vector to hold the values
        std::vector<double> values;
        double value;

        // parse the line
        while ((pos = line.find(delimiter)) != std::string::npos) {
            std::string token = line.substr(0, pos);
            value = std::stod(token);
            values.push_back(value);
            line.erase(0, pos + delimiter.length());
        }

        // add the last value
        if (!line.empty()) {
            value = std::stod(line);
            values.push_back(value);
        }

        // create a DataPoint object and add it to the points vector
        if (!values.empty()) {
            double* data = new double[values.size()];
            std::copy(values.begin(), values.end(), data);
            DataPoint pt(data, values.size());
            points.push_back(pt);
            delete[] data; // Free the allocated memory for data
        }
    }

    file.close();
    return points.size();
}
