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

struct Rect {
    Rect() {}

    double min[2];
    double max[2];

    Rect(double a_minX, double a_minY,  
         double a_maxX, double a_maxY) {
      min[0] = a_minX;
      min[1] = a_minY;
  
      max[0] = a_maxX;
      max[1] = a_maxY;
    }
};

DBSCAN::DBSCAN(vector<DataPoint>& points, int minPts, double eps) {
    // Initialize parameters
    this->minPts = minPts;
    this->eps = eps;
    this->cluster_id = 0;
    this->data_size = points.size();
    this->data = new vector<DataPoint>(points);
    this->labels = new int[this->data_size](); // Allocate memory for labels
    
    // Dynamically determine the dimension from the first DataPoint
    if (!points.empty()) {
        this->dim = points[0].get_dim();
    } else {
        this->dim = 0;
    }

    // Allocate the dataset matrix for FLANN
    this->dataset = flann::Matrix<double>(new double[this->data_size * this->dim], this->data_size, this->dim);

    // Populate the dataset matrix with the data from the points
    for (size_t i = 0; i < this->data_size; i++) {
        for (int j = 0; j < this->dim; j++) {
            this->dataset[i][j] = (*this->data)[i][j];
        }
    }

    // Build the FLANN k-d tree index
    this->index = new flann::Index<flann::L2<double>>(this->dataset, flann::KDTreeIndexParams(4));
    this->index->buildIndex();
}

DBSCAN::~DBSCAN() {
    delete this->data; // Free the allocated memory for data points
    delete[] labels; // Free the allocated memory for labels
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
vector<int> DBSCAN::regionQuery(int point, vector<DataPoint>& points) {
    vector<int> neighbors;

    // Prepare the query point
    flann::Matrix<double> query(new double[this->dim], 1, this->dim);
    for (int i = 0; i < this->dim; i++) {
        query[0][i] = (*this->data)[point][i];
    }

    // Perform a radius search
    std::vector<std::vector<int>> indices;
    std::vector<std::vector<double>> dists;
    this->index->radiusSearch(query, indices, dists, eps * eps, flann::SearchParams(128));

    // Add the neighbors to the result
    if (!indices.empty()) {
        neighbors = indices[0]; // Get the neighbors for the first query point
    }

    delete[] query.ptr(); // Free the query matrix memory
    return neighbors;
}

void DBSCAN::run() {
    // Iterate through points
    for (int i = 0; i < this->data_size; i++) {
        if (this->labels[i] != 0) {
            continue; // Already processed
        }
        // Search for neighbors within eps distance
        vector<int> neighbors = regionQuery(i, *this->data);
        if (neighbors.size() < this->minPts) {
            printf("Point %d is noise\n", i);
            // Mark as noise
            this->labels[i] = NOISE; // Mark as noise
            continue;
        }

        this->cluster_id++; // New cluster found

        // increment the cluster id for the current point
        this->labels[i] = this->cluster_id;

        // Expand the cluster
        while (!neighbors.empty()) {
            // get the neighbor index
            int neighbor = neighbors.front();

            if (neighbor == i) {
                continue; // Skip the point itself
            }

            // Change noise to cluster id
            if (this->labels[neighbor] == NOISE) {
                this->labels[neighbor] = this->cluster_id; 
                continue;
            }

            if (this->labels[neighbor] != 0) {
                continue; // Already processed
            }

            this->labels[neighbor] = this->cluster_id; // Assign cluster id
            
            // Search for neighbors of the neighbor
            vector<int> newNeighbors = regionQuery(neighbor, *this->data);
            if (newNeighbors.size() >= this->minPts) {

                // Merge the new neighbors into the existing neighbors
                neighbors.insert(neighbors.end(), newNeighbors.begin(), newNeighbors.end());
            }

        }
    }
}

void DBSCAN::results() {
    printf("Number of clusters: %d\n", this->cluster_id);
    int noises = 0;
    for (size_t i = 0; i < this->data_size; i++) {
      if (this->labels[i] == NOISE) {
        noises++;
      }
    }
  
    printf("Noises: %d\n", noises);
}
