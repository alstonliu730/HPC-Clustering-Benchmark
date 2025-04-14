#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>
#include "datapoint.h"
#include <cstddef> // for size_t

size_t import_data(char* filename, std::vector<DataPoint>& points);

#endif // _UTILS_H_
