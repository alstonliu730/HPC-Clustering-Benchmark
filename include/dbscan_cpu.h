#ifndef _DBSCAN_CPU_H_
#define _DBSCAN_CPU_H_

#include <vector>
#include <optional>
#include "datapoint.h"
#include "flann/flann.hpp"

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#define NOISE -1

using namespace flann;
class DBSCAN {
    protected:
        int minPts; // Minimum number of points in a neighborhood to form a dense region
        double eps; // Maximum distance between two points to be considered neighbors

        size_t cluster_id; // Current cluster ID
        std::vector<size_t>* labels; // Array to store cluster labels for each point
        int dim; // Dimension of the data points

        std::vector<DataPoint>* data; // pointer to the data points
        size_t data_size; // Size of the data array

        Matrix<float> dataset; // Matrix to store the dataset for FLANN
        KDTreeSingleIndex<L2_Simple<float>>* index; // FLANN index for nearest neighbor search

    private:
        // Function to find all points within eps distance from the given point
        std::vector<size_t> regionQuery(size_t p_idx, const DataPoint& point, const int max_nn = 10);
        #ifdef USE_MPI
        void serializeKDTree(flann::KDTreeSingleIndex<flann::L2<float>>& index, std::vector<char>& buffer);
        void broadcastKDTree(flann::KDTreeSingleIndex<flann::L2<float>>& index, int root, MPI_Comm comm);
        #endif // USE_MPI
    public:
        DBSCAN(const std::vector<DataPoint>& points, int minpts, double eps);
        ~DBSCAN();
        void run();
        void result(std::vector<size_t>& labels);
        
        #ifdef USE_MPI
        void mpi_run();
        #endif // USE_MPI
};

#endif // _DBSCAN_CPU_H_
