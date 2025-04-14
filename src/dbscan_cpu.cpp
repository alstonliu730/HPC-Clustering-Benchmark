#include <iostream>
#include <vector>
#include <math.h>
#include <omp.h>
#include <functional>
#include "dbscan_cpu.h"
#include "datapoint.h"
#include "rtree.h"
#include "utils.h"

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
    this->data = new std::vector<DataPoint>(points);
    this->labels = new int[this->data_size]; // Allocate memory for labels
    
    // Build the RTree for spatial indexing
    for (size_t i = 0; i < this->data_size; ++i) {
        Rect rect(points[i][0] - eps, points[i][1] - eps,
                  points[i][0] + eps, points[i][1] + eps);
        tree.Insert(rect.min, rect.max, i);
    }
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
        std::cerr << "Error: Dimensions of the points do not match." << std::endl;
        return -1.0; // Return an error value
    }

    double dist = 0.0;
    int i;
    for(i = 0; i < a.get_dim(); i++) {
        double diff = a[i] - b[i];
        dist += diff * diff; // Squared distance
    }

    return sqrt(dist); // Return the Euclidean distance
}

/**
 * @brief Find all points within eps distance from the given point
 * @return A vector of indices of the points within eps distance from the given point
 */
vector<int> DBSCAN::regionQuery(int point, vector<DataPoint>& points) {
    // Create a vector to hold the neighbors
    vector<int> neighbors, searchNeighbors;

    // Create a search rectangle around the point
    DataPoint p = points[point];    
    Rect searchRect = Rect(p[0] - eps, p[1] - eps,
                           p[0] + eps, p[1] + eps);

    // Create a callback function to add the found points to the neighbors                      
    std::function<bool (const int)> searchBoxCallback = [&](const int id) {
        searchNeighbors.push_back(id);
        return true;
    };

    // Search the RTree for points within the search rectangle
    this->tree.Search(searchRect.min, searchRect.max, searchBoxCallback);
    
    // Iterate through the neighbors and check if they are within eps distance
    #pragma omp parallel 
    {   
        int chunk_size = (searchNeighbors.size() / omp_get_num_threads()) + 1;  1;
        #pragma omp for schedule(static, chunk_size)
        for (int i = 0; i < searchNeighbors.size(); i++) {
            DataPoint neighbor = points[searchNeighbors[i]];
            if (getDist(p, neighbor) <= eps) {
                #pragma omp critical
                neighbors.push_back(searchNeighbors[i]);
            }
        }
    }

    // return the neighbors
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

            if (this->labels[neighbor] == NOISE) {
                this->labels[neighbor] = this->cluster_id; // Change noise to cluster id
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