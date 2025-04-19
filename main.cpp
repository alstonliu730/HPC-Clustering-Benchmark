#include "utils.h"
#include "datapoint.h"
#include "dbscan_cpu.h"
#include "flann/flann.hpp"
#include <iostream>
#include <string>
#include <chrono>
#include <vector>

using namespace std;

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        return 1;
    }

    char* input_file = argv[1];
    vector<DataPoint> points;

    size_t data_size = import_data(input_file, points);
    if (data_size == 0) {
        cerr << "Error: No data points found in the file." << endl;
        return 1;
    }

    cout << "Number of data points imported: " << data_size << endl;
    // Start the timer
    auto start = chrono::high_resolution_clock::now();

    // Trying high eps distance
    DBSCAN dbscan(points, 5, 1.0f); // Example parameters: minPts = 5, eps = 1
    dbscan.run(); // Run the DBSCAN algorithm
    
    // Benchmarking the DBSCAN algorithm
    auto end = chrono::high_resolution_clock::now();
    auto ms = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    cout << "DBSCAN completed in " << ms << " ms" << endl;

    vector<size_t> labels;
    dbscan.result(labels); // Get the labels for each point
    
    // Export the results to a file
    string output_file = "output.csv";
    size_t exported_size = export_data(output_file.c_str(), points, labels);
    if (exported_size == 0) {
        cerr << "Error: No data points exported to the file." << endl;
        return 1;
    }

    dbscan.~DBSCAN(); // Clean up the DBSCAN object
    return 0;
}
