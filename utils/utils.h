#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>
#include "datapoint.h"
#include <cstddef> // for size_t

size_t import_data(char* filename, std::vector<DataPoint>& points);
size_t export_data(const char* output_file, const std::vector<DataPoint>& points,
    const std::vector<size_t>& labels);

#endif // _UTILS_H_
