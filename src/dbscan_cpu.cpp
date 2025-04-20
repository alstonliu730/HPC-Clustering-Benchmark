#include <iostream>
#include <vector>
#include <math.h>
#include <omp.h>
#include <sstream>
#include <vector>

#include "dbscan_cpu.h"
#include "datapoint.h"
#include "utils.h"
#include "flann/flann.hpp"

#ifdef USE_MPI
#include <mpi.h>

#define ROOT 0

void DBSCAN::serializeKDTree(flann::KDTreeSingleIndex<flann::L2<float>>& index, std::vector<char>& buffer) {
    std::ostringstream oss(std::ios::binary);
    flann::serialization::SaveArchive sa(oss);
    sa & index; // Serialize the k-d tree into the stream
    std::string serialized_data = oss.str();
    buffer.assign(serialized_data.begin(), serialized_data.end()); // Copy to buffer
}

void DBSCAN::broadcastKDTree(flann::KDTreeSingleIndex<flann::L2<float>>& index, MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::vector<char> buffer;

    if (rank == ROOT) {
        // Serialize the KDTree on the root process
        serializeKDTree(index, buffer);
    }

    // Broadcast the size of the serialized data
    int buffer_size = buffer.size();
    MPI_Bcast(&buffer_size, 1, MPI_INT, ROOT, comm);

    // Resize the buffer on non-root processes
    if (rank != ROOT) {
        buffer.resize(buffer_size);
    }

    // Broadcast the serialized data
    MPI_Bcast(buffer.data(), buffer_size, MPI_CHAR, ROOT, comm);

    // Deserialize the KDTree on non-root processes
    if (rank != ROOT) {
        std::istringstream iss(std::string(buffer.begin(), buffer.end()), std::ios::binary);
        flann::serialization::LoadArchive la(iss);
        la & index; // Deserialize the k-d tree
    }
}

void DBSCAN::mpi_run() {
    // Initialize MPI
    MPI_Init(nullptr, nullptr);

    // Get the number of processes and the rank of this process
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Calculate the local size for each process
    std::vector<int> sendcounts(world_size);
    std::vector<int> displs(world_size);
    std::vector<float> sendbuffer(this->data_size * this->dim); // flatten the data points
    if (world_rank == 0) {
        // Prepare the send buffer for the root process
        for (size_t i = 0; i < this->data_size; i++) {
            for (int j = 0; j < this->dim; j++) {
                sendbuffer[i * this->dim + j] = (*this->data)[i][j];
            }
        }

        // Calculate the send counts and displacements for each process
        printf("Distributing data points among %d processes...\n", world_size);
        int total_size = this->data_size * this->dim;
        int offset = 0;
        for (int i = 0; i < world_size; i++) {
            sendcounts[i] = total_size / world_size 
                + (i < total_size % world_size ? 1 : 0);
            displs[i] = offset;
            offset += sendcounts[i];
        }
    }
    // Broadcast the dimension of the data points to all processes
    MPI_Bcast(&this->dim, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&this->data_size, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&this->minPts, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&this->eps, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    // Tell each process how many items it will receive
    int local_count; // Number of floating point numbers to receive
    MPI_Scatter(sendcounts.data(), 1, MPI_INT, &local_count, 1, MPI_INT, ROOT, MPI_COMM_WORLD);   

    vector<float> local_data(local_count); // Allocate memory for local data points (x0, y0, z0, ...)
    MPI_Scatterv(sendbuffer.data(), sendcounts.data(), displs.data(), MPI_FLOAT, 
        local_data.data(), local_count, MPI_FLOAT, ROOT, MPI_COMM_WORLD);

    // Broadcast the KDTree from the root process to all other processes
    broadcastKDTree(*this->index, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD); // Synchronize all processes

    // DEBUGGING 
    if (world_rank == 1) {
        // Check the received data
        printf("Process %d received %d data points:\n", world_rank, local_count / this->dim);
        
        // Check if the index is valid
        this->index->
    }
    // Create a local vector of DataPoint objects for each process
    size_t local_size = local_count / this->dim; // Number of data points for this process
    vector<DataPoint> local_points(local_size);

    // Fill the local points with the received data
    for (size_t i = 0; i < local_size; i++) {
        float* data_ptr = new float[this->dim]; // Allocate memory for the data array
        for (int j = 0; j < this->dim; j++) {
            data_ptr[j] = local_data[i * this->dim + j];
        }
        DataPoint point(data_ptr, this->dim); // Create a DataPoint object
        local_points[i] = point; // Assign the DataPoint object to the local points vector
        delete[] data_ptr; // Free the allocated memory for the data array
    }

    // Create a local labels vector for each process
    vector<int> local_labels(local_size, -1); // Initialize labels to 0
    for (size_t i = 0; i < local_size; i++) {
        local_labels[i] = 0; // Initialize labels to 0
    }
}
#endif // USE_MPI

using namespace std;

/**
 * @brief Constructor for DBSCAN class
 * @param points Vector of DataPoint objects
 * @param minPts Minimum number of points in a neighborhood to form a dense region
 * @param eps Maximum distance between two points to be considered neighbors
 */
DBSCAN::DBSCAN(const vector<DataPoint>& points, int minPts, double eps) {
    // Initialize parameters
    printf("----- Initializing DBSCAN -----\n");
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
    this->labels = new vector<size_t>(this->data_size); // Allocate memory for labels
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
    this->index = new flann::KDTreeSingleIndex<flann::L2_Simple<float>>(this->dataset, flann::KDTreeSingleIndexParams(10));
    this->index->buildIndex();
}

/**
 * @brief Destructor for DBSCAN class
 */
DBSCAN::~DBSCAN() {
    printf("----- Cleaning up -----\n");
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
vector<size_t> DBSCAN::regionQuery(const size_t p_idx, const DataPoint& point, const int max_nn = 10) {
    // printf("Finding neighbors for point %ld\n", point);
    vector<size_t> neighbors;

    // Prepare the query point
    flann::Matrix<float> query(new float[this->dim], 1, this->dim);
    for (int i = 0; i < this->dim; i++) {
        query[0][i] = point[i];
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

    for (size_t i = 0; (i < num_found && i < max_nn); i++) {
        size_t idx = indices[0][i];
        if (idx != p_idx) { // Skip the point itself
            neighbors.push_back(idx); // Add the neighbor index to the vector
            // printf("Neighbor %ld: %ld\n", i, indices[i][0]);
        }
    }

    delete[] query.ptr(); // Free the query matrix memory
    delete[] indices.ptr(); // Free the indices matrix memory
    delete[] dists.ptr(); // Free the distances matrix memory
    return neighbors;
}

/**
 * @brief Run the DBSCAN algorithm on the class data
 */
void DBSCAN::run() {
    int max_nn = this->minPts * 2; // Maximum number of neighbors to search for
    printf("------ Running DBSCAN with max_nn(%d) & minPts(%d) ------\n", max_nn, this->minPts);

    // Iterate through points
    for (size_t i = 0; i < this->data_size; i++) {
        // printf("Processing point %ld\n", i);
        if ((*this->labels)[i] != 0) continue;

        vector<size_t> neighbors = regionQuery(i, (*this->data)[i]);
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
            vector<size_t> newNeighbors = regionQuery(neighbor, (*this->data)[neighbor]);
            if (newNeighbors.size() >= this->minPts) {
                // Merge the new neighbors into the existing neighbors
                neighbors.insert(neighbors.end(), newNeighbors.begin(), newNeighbors.end());
            }
        }
    }
}

/**
 * @brief Get the result of the clustering
 * @param res Vector to store the cluster labels for each point
 * @return void
 */
void DBSCAN::result(std::vector<size_t>& res) {
    printf("------ DBSCAN Result ------\n");
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
