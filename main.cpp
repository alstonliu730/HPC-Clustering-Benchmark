#include "utils.h"
#include "datapoint.h"
#include "dbscan_cpu.h"
#include "rtree.h"
#include "flann/flann.hpp"
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1;
    }

    char* input_file = argv[1];
    std::vector<DataPoint> points;

    size_t data_size = import_data(input_file, points);
    if (data_size == 0) {
        std::cerr << "Error: No data points found in the file." << std::endl;
        return 1;
    }

    std::cout << "Number of data points imported: " << data_size << std::endl;

    DBSCAN dbscan(points, 5, 0.5); // Example parameters: minPts = 5, eps = 0.5
    dbscan.run(); // Run the DBSCAN algorithm
    dbscan.results(); // Output the results
}
