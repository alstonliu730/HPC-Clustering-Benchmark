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
    this->data = vector<DataPoint>(points.begin(), points.end()); // Copy the data points
    this->labels = vector<std::atomic<size_t>>(this->data_size); // Allocate memory for labels
    // this->labels->assign(this->data_size, 0); // Initialize labels to 0
    this->visited = vector<bool>(this->data_size, false); // Initialize visited array to false
    
    // Allocate the dataset matrix for FLANN
    this->dataset = flann::Matrix<float>(new float[this->data_size * this->dim], this->data_size, this->dim);

    // Populate the dataset matrix with the data from the points
    printf("Populating dataset matrix...\n");
    omp_set_num_threads(4); // Set number of threads before the parallel block
    omp_set_nested(0);
    #pragma omp parallel 
    {
    int chunk_size = this->data_size / omp_get_num_threads(); // Calculate chunk size for parallel processing
    #pragma omp for schedule(static, chunk_size)
    for (size_t i = 0; i < this->data_size; i++) {
        for (int j = 0; j < this->dim; j++) {
            this->dataset[i][j] = this->data[i][j];
        }
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
    this->labels.clear(); // Clear the labels vector

    printf("Deleting data...\n");
    this->data.clear(); // Clear the data vector

    printf("Deleting FLANN index...\n");
    this->index->~KDTreeSingleIndex(); // Free the FLANN index

    printf("Deleting dataset...\n");
    delete[] this->dataset.ptr(); // Free the dataset matrix memory

    printf("Clearing visited...\n");
    this->visited.clear();
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
    int num_found = this->index->radiusSearch(query, indices, dists, (this->eps * this->eps), flann::SearchParams(64, 0.f, false));
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
    const vector<DataPoint> data = this->data; // Copy the data points
    atomic<size_t> next_cluster_id(1); // Atomic variable to ensure thread safety for cluster ID

    // Iterate through points
    #pragma omp parallel 
    {
        // Initialize OpenMP parameters
        size_t nThreads = omp_get_num_threads(); // Number of threads available
        // printf("Number of threads: %ld\n", nThreads);
        size_t chunk_size = (nPoints / nThreads) + 1; // Calculate chunk size for parallel processing

	    #pragma omp single 
        {
            printf("DBSCAN threads: %ld\n", nThreads);
            printf("Chunk size: %ld\n", chunk_size);
        }
        
        #pragma omp parallel for schedule(static, chunk_size)
        for (size_t i = 0; i < nPoints; i++) {
            if (this->visited[i]) continue;
            this->visited[i] = true; // Mark the point as visited

            vector<size_t> neighbors = regionQuery(i, data);
            if (neighbors.size() < static_cast<size_t>(minPts)) {
                this->labels[i].store(NOISE, memory_order_release);
                continue;
            }
            
            size_t local_cluster_id;
            local_cluster_id = next_cluster_id.fetch_add(1, memory_order_relaxed); // Increment the cluster ID atomically
    
            this->labels[i].store(local_cluster_id, memory_order_release); // Assign the cluster ID to the current point
    
            vector<size_t> stack(neighbors.begin(), neighbors.end());
            while (!stack.empty()) {
                size_t neighbor = stack.back();
                stack.pop_back();
    
                if (neighbor == i) continue;
                if (this->visited[neighbor]) continue;

                this->visited[neighbor] = true; // Mark the point as visited

                size_t prev = this->labels[neighbor].load(memory_order_acquire);
                if (prev == NOISE) {
                    this->labels[neighbor].store(local_cluster_id, memory_order_release);
                    continue;
                }
    
                bool updated = this->labels[neighbor].compare_exchange_strong(prev, local_cluster_id, memory_order_relaxed);
                if (updated) {
                    vector<size_t> new_neighbors = regionQuery(neighbor, data);
                    if (new_neighbors.size() >= static_cast<size_t>(minPts)) {
                        copy(new_neighbors.begin(), new_neighbors.end(), back_inserter(stack)); // Add new neighbors to the list
                    }
                }
            }

            // clear the stack for the next iteration
            stack.clear();
            neighbors.clear();
        }

        #pragma omp barrier // Ensure all threads have completed before updating the cluster ID
        #pragma omp single
        {
            this->cluster_id = next_cluster_id.load(memory_order_relaxed) - 1; // Update the cluster ID
        }
    }
}

void DBSCAN::result(std::vector<size_t>& res) {
    printf("Number of clusters: %ld\n", this->cluster_id);
    size_t noises = 0;
    for (size_t i = 0; i < this->data_size; i++) {
      if (this->labels[i].load() == NOISE) {
        noises++;
      }
    }
  
    printf("Noises: %ld\n", noises);
    res.assign(this->labels.begin(), this->labels.end());
}
