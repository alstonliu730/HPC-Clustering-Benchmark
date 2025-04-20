#include <iostream>
#include <vector>
#include <math.h>
#include <omp.h>
#include <functional>
#include <algorithm>
#include <atomic>

#include "dbscan_cpu.h"
#include "datapoint.h"
#include "rtree.h"
#include "utils.h"
#include "flann/flann.hpp"

using namespace std;

DBSCAN::DBSCAN(const vector<DataPoint>& points, int minPts, double eps) {
    // Initialize parameters
    this->minPts = minPts;
    this->eps = eps;
    this->cluster_id = 0;
    this->data_size = points.size();

    // Dynamically determine the dimension from the first DataPoint
    if (!points.empty()) {
        this->dim = points[0].get_dim();
    } else {
        this->dim = 0;
    }
    this->data = new vector<DataPoint>(points.begin(), points.end()); // Copy the data points
    this->labels = new vector<std::atomic<size_t>>(this->data_size); // Allocate memory for labels
    this->labels->assign(this->data_size, 0); // Initialize labels to 0

    // Allocate the dataset matrix for FLANN
    this->dataset = flann::Matrix<float>(new float[this->data_size * this->dim], this->data_size, this->dim);

    // Populate the dataset matrix with the data from the points
    printf("Populating dataset matrix...\n");
    int chunk_size = this->data_size / omp_get_num_threads(); // Calculate chunk size for parallel processing
    #pragma omp parallel for schedule(static, chunk_size)
    for (size_t i = 0; i < this->data_size; i++) {
        for (int j = 0; j < this->dim; j++) {
            this->dataset[i][j] = (*this->data)[i][j];
        }
    }

    // Build the FLANN k-d tree index
    printf("Building FLANN index...\n");
    this->index = new flann::KDTreeSingleIndex<flann::L2_Simple<float>>(this->dataset, flann::KDTreeSingleIndexParams(16));
    this->index->buildIndex();
}

DBSCAN::~DBSCAN() {
    printf("Cleaning up...\n");
    printf("Deleting labels...\n");
    this->labels->clear(); // Clear the labels vector
    delete this->labels; // Free the labels array
    this->labels = nullptr; // Set to nullptr to avoid dangling pointer

    printf("Deleting data...\n");
    this->data->clear(); // Clear the data vector
    delete this->data; // Free the data vector
    this->data = nullptr; // Set to nullptr to avoid dangling pointer

    printf("Deleting FLANN index...\n");
    this->index->~KDTreeSingleIndex(); // Free the FLANN index

    printf("Deleting dataset...\n");
    delete[] this->dataset.ptr(); // Free the dataset matrix memory
}

/**
 * @brief Find all points within eps distance from the given point
 * @return A vector of indices of the points within eps distance from the given point
 */
vector<size_t> DBSCAN::regionQuery(size_t point, const vector<DataPoint>& points) {
    // printf("Finding neighbors for point %ld\n", point);
    vector<size_t> neighbors;
    int max_nn = this->minPts * 5; // Maximum number of neighbors to search for

    // Prepare the query point
    flann::Matrix<float> query(new float[this->dim], 1, this->dim);
    for (int i = 0; i < this->dim; i++) {
        query[0][i] = points[point][i];
    }

    // Perform a radius search
    // printf("Searching for neighbors...\n");
    flann::Matrix<size_t> indices(new size_t[max_nn], 1, max_nn); // Allocate memory for indices
    flann::Matrix<float> dists(new float[max_nn], 1, max_nn); // Allocate memory for distances

    // Perform the radius search using FLANN
    int num_found = this->index->radiusSearch(query, indices, dists, (this->eps * this->eps), flann::SearchParams(32));
    if (num_found == -1) {
        printf("Error: No neighbors found.\n");
        return neighbors; // Return empty vector if no neighbors found
    }

    //printf("Neighbors found: %d\n", num_found); // Print the number of neighbors found

    for (size_t i = 0; i < num_found; i++) {
        size_t idx = indices[0][i];
        if (idx != point && neighbors.size() < max_nn) { // Skip the point itself
            neighbors.push_back(idx); // Add the neighbor index to the vector
            // printf("Neighbor %ld: %ld\n", i, indices[i][0]);
        }
    }

    delete[] query.ptr(); // Free the query matrix memory
    delete[] indices.ptr(); // Free the indices matrix memory
    delete[] dists.ptr(); // Free the distances matrix memory
    return neighbors;
}

void DBSCAN::run() {
    printf("Running DBSCAN...\n");
    // Initialize parameters for data locality
    int minPts = this->minPts; // Minimum number of points in a neighborhood to form a dense region
    size_t nPoints = this->data_size; // Number of data points
    size_t next_cluster_id = 1;
    const vector<DataPoint> data = *this->data; // Copy the data points
    vector<std::atomic<size_t>> labels = *this->labels; // Copy the labels

    // Initialize OpenMP parameters
    size_t nThreads = omp_get_max_threads() / 2; // Number of threads available
    printf("Number of threads: %ld\n", nThreads);
    size_t chunk_size = (nPoints / nThreads) + 1; // Calculate chunk size for parallel processing

    // Iterate through points
    #pragma omp parallel for schedule(static, chunk_size) shared(next_cluster_id)
    for (size_t i = 0; i < nPoints; i++) {
        if (labels[i].load() != 0) continue;

        vector<size_t> neighbors = regionQuery(i, data);
        if (neighbors.size() < minPts) {
            labels[i].store(NOISE);
            continue;
        }

        size_t local_cluster_id;
        #pragma omp critical
        {
            local_cluster_id = next_cluster_id++;
        }

        labels[i].store(local_cluster_id);

        vector<size_t> stack(neighbors.begin(), neighbors.end());
        while (!stack.empty()) {
            size_t neighbor = stack.back();
            stack.pop_back();

            if (neighbor == i) continue;

            size_t prev = labels[neighbor].load();
            if (prev == NOISE) {
                labels[neighbor].store(local_cluster_id);
                continue;
            }

            if (prev != 0) continue;

            bool updated = labels[neighbor].compare_exchange_strong(prev, local_cluster_id);
            if (updated) {
                vector<size_t> new_neighbors = regionQuery(neighbor, data);
                if (new_neighbors.size() >= minPts) {
                    stack.insert(stack.end(), new_neighbors.begin(), new_neighbors.end());
                }
            }
        }
    }

    // Update the labels vector
    this->labels->assign(labels.begin(), labels.end());
    this->cluster_id = next_cluster_id - 1; // Update the cluster ID
}

void DBSCAN::result(std::vector<size_t>& res) {
    printf("Number of clusters: %ld\n", this->cluster_id);
    size_t noises = 0;
    for (size_t i = 0; i < this->data_size; i++) {
      if ((*this->labels)[i].load() == NOISE) {
        noises++;
      }
    }
  
    printf("Noises: %ld\n", noises);
    res.assign(this->labels->begin(), this->labels->end());
}
