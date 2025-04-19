#include "utils.h"
#include "datapoint.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <iomanip>

size_t import_data(char* input_file, std::vector<DataPoint>& points) {
    std::ifstream file(input_file);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << input_file << std::endl;
        return 0;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::vector<float> values;
        std::string token;

        while (std::getline(ss, token, ',')) {
            if (!token.empty()) {
                values.push_back(std::stod(token));
            }
        }

        if (!values.empty()) {
            float* data = new float[values.size()];
            std::copy(values.begin(), values.end(), data);
            DataPoint pt(data, values.size());
            points.push_back(pt);
            delete[] data;
        }
    }

    file.close();
    return points.size();
}

size_t export_data(const char* output_file, const std::vector<DataPoint>& points,
    const std::vector<size_t>& labels) {
    assert(points.size() == labels.size());
    std::ofstream file(output_file);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << output_file << std::endl;
        return 0;
    }

    file << std::fixed << std::setprecision(6);

    for (int i = 0; i < points.size(); i++) {
        DataPoint point = points[i];
        for (int j = 0; j < point.get_dim(); j++) {
            file << point[j] << ","; // Write the data point
        }
	size_t s = labels[i];
        if (s == size_t(-1)) {
            s = 0;
	}
        file << s; // Append the label  
        file << "\n";
    }

    file.close();
    return points.size();
}
