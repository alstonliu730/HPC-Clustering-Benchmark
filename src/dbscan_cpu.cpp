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
    printf("Populating dataset matrix...\n");
    #pragma omp parallel for schedule(static, (data_size/omp_get_num_threads())) // Parallelize the loop for performance
    for (size_t i = 0; i < this->data_size; i++) {
        for (int j = 0; j < this->dim; j++) {
            this->dataset[i][j] = (*this->data)[i][j];
        }
    }

    // Build the FLANN k-d tree index
    printf("Building FLANN index...\n");
    this->index = new flann::Index<flann::L2<double>>(this->dataset, flann::KDTreeIndexParams(4));
    this->index->buildIndex();
}

DBSCAN::~DBSCAN() {
    delete[] this->labels; // Free the labels array
    delete this->data; // Free the data vector
    this->index->~Index(); // Free the FLANN index
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
vector<int> DBSCAN::regionQuery(int point, vector<DataPoint>& points) {
    printf("Finding neighbors for point %d\n", point);
    vector<int> neighbors;

    // Prepare the query point
    flann::Matrix<double> query(new double[this->dim], 1, this->dim);
    for (int i = 0; i < this->dim; i++) {
        query[0][i] = points[point][i];
    }

    // Perform a radius search
    printf("Searching for neighbors...\n");
    flann::Matrix<int> indices; // Allocate memory for indices
    flann::Matrix<double> dists; // Allocate memory for distances

    // Perform the radius search using FLANN
    int num_found = this->index->radiusSearch(query, indices, dists, (this->eps * this->eps), flann::SearchParams(32));
    if (num_found == -1) {
        printf("Error: No neighbors found.\n");
        return neighbors; // Return empty vector if no neighbors found
    }

    printf("Neighbors found: %d\n", num_found); // Print the number of neighbors found
    
    // Add the neighbors to the result
    int* indices_res = indices.ptr();
    double* dists_res = dists.ptr();
    std::cout << "Neighbors within radius " << this->eps << ":\n";
    for (int i = 0; i < indices.rows; i++) {
        for(int j = 0; j < indices.cols; j++) {
            std::cout << " - Point index: " << indices_res[i*this->dim + j] << '\n';
            neighbors.push_back(indices_res[i*this->dim + j]); // Add the neighbor index to the result
        }
    }

    delete[] query.ptr(); // Free the query matrix memory
    delete[] indices.ptr(); // Free the indices matrix memory
    delete[] dists.ptr(); // Free the distances matrix memory
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
            this->labels[i] = NOISE; 
            continue;
        }

        this->cluster_id++; // New cluster found

        // increment the cluster id for the current point
        this->labels[i] = this->cluster_id;

        // Expand the cluster
        for(int neighbor : neighbors) {
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
