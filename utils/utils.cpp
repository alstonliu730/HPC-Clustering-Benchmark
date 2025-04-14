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
        std::stringstream ss(line);
        std::vector<double> values;
        std::string token;

        while (std::getline(ss, token, ',')) {
            if (!token.empty()) {
                values.push_back(std::stod(token));
            }
        }

        if (!values.empty()) {
            double* data = new double[values.size()];
            std::copy(values.begin(), values.end(), data);
            DataPoint pt(data, values.size());
            points.push_back(pt);
            delete[] data;
        }
    }

    file.close();
    return points.size();
}
