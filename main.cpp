#include "utils.h"
#include "datapoint.h"
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

    std::cout << "Successfully loaded " << data_size << " data points." << std::endl;

    // Print a few points to verify
    for (int i = 0; i < std::min(data_size, size_t(5)); ++i) {
        std::cout << "Point " << i << ": ";
        for (int j = 0; j < points[i].get_dim(); ++j)
            std::cout << points[i][j] << " ";
        std::cout << "\n";
    }

    return 0;
}
