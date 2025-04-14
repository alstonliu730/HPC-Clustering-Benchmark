#include "dbscan_cpu.h"
#include "utils.h"
#include "datapoint.h"
#include <iostream>
#include <string.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <minPts> <eps>" << std::endl;
        return 1;
    }

    char* input_file = argv[1];
    int minPts = std::stoi(argv[2]);
    double eps = std::stod(argv[3]);

    std::vector<DataPoint> points;
    size_t data_size = import_data(input_file, points);
    if (data_size == 0) {
        std::cerr << "Error: No data points found in the file." << std::endl;
        return 1;
    }
    std::cout << "Data size: " << data_size << std::endl;
    return 0;
}
