#include <iostream>
#include <vector>
#include <math.h>
#include <omp.h>
#include <functional>
#include "dbscan_cpu.h"
#include "datapoint.h"
#include "rtree.h"
#include "utils.h"
#include "flann/flann.hpp"

using namespace std;

DBSCAN::DBSCAN(vector<DataPoint>& points, int minPts, double eps) {
    // Initialize parameters
    this->minPts = minPts;
    this->eps = eps;
    this->cluster_id = 0;
    this->data_size = points.size();
    this->data = new vector<DataPoint>(points);
    this->labels = new vector<size_t>(this->data_size); // Allocate memory for labels
    
    // Dynamically determine the dimension from the first DataPoint
    if (!points.empty()) {
        this->dim = points[0].get_dim();
    } else {
        this->dim = 0;
    }

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
    this->index = new flann::KDTreeSingleIndex<flann::L2_Simple<float>>(this->dataset, flann::KDTreeSingleIndexParams(10));
    this->index->buildIndex();
}

DBSCAN::~DBSCAN() {
    printf("Cleaning up...\n");
    printf("Deleting labels...\n");
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
 * @brief Return the euclidean distance between two datapoints.
 */
double DBSCAN::getDist(DataPoint &a, DataPoint &b) {
    if (a.get_dim() != b.get_dim()) {
        cerr << "Error: Dimensions of the points do not match." << std::endl;
        return -1.0; // Return an error value
    }

    double dist = 0.0;
    int i;
    for(i = 0; i < a.get_dim(); i++) {
        double diff = a[i] - b[i];
        dist += (diff * diff); // Squared distance
    }

    return sqrt(dist); // Return the Euclidean distance
}

/**
 * @brief Find all points within eps distance from the given point
 * @return A vector of indices of the points within eps distance from the given point
 */
vector<size_t> DBSCAN::regionQuery(size_t point, vector<DataPoint>& points) {
    // printf("Finding neighbors for point %ld\n", point);
    vector<size_t> neighbors;
    int max_nn = 10;
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

    printf("Neighbors found: %d\n", num_found); // Print the number of neighbors found
    
    for (size_t i = 0; i < num_found; i++) {
        size_t idx = indices[0][i];
        if (idx != point && i < max_nn) { // Skip the point itself
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
    // Iterate through points
    for (size_t i = 0; i < this->data_size; i++) {
        // printf("Processing point %ld\n", i);
        if ((*this->labels)[i] != 0) continue;

        vector<size_t> neighbors = regionQuery(i, *this->data);
        if (neighbors.size() < this->minPts) {
            (*this->labels)[i] = NOISE;
            continue;
        }
        // printf("Point %ld is a core point\n", i);

        // Assign a new cluster id
        this->cluster_id++;
        (*this->labels)[i] = this->cluster_id;

        // Expand the cluster
        // printf("Expanding cluster %ld\n", this->cluster_id);
        while (!neighbors.empty()) {
            size_t neighbor = neighbors.back();
            neighbors.pop_back();
            // printf("Processing neighbor %ld\n", neighbor);
            if (neighbor == i) {
                continue; // Skip the point itself
            }

            // Change noise to cluster id
            if ((*this->labels)[neighbor] == NOISE) {
                // printf("Point %ld is no longer noise\n", neighbor);
                (*this->labels)[neighbor] = this->cluster_id; 
                continue;
            }

            if ((*this->labels)[neighbor] != 0) {
                continue; // Already processed
            }

            (*this->labels)[neighbor] = this->cluster_id; // Assign cluster id
            
            // Search for neighbors of the neighbor
            vector<size_t> newNeighbors = regionQuery(neighbor, *this->data);
            if (newNeighbors.size() >= this->minPts) {
                // Merge the new neighbors into the existing neighbors
                neighbors.insert(neighbors.end(), newNeighbors.begin(), newNeighbors.end());
            }
        }
    }
}

void DBSCAN::result(std::vector<size_t>& res) {
    printf("Number of clusters: %ld\n", this->cluster_id);
    size_t noises = 0;
    for (size_t i = 0; i < this->data_size; i++) {
      if ((*this->labels)[i] == NOISE) {
        noises++;
      }
    }
  
    printf("Noises: %ld\n", noises);
    res.assign(this->labels->begin(), this->labels->end());
}
