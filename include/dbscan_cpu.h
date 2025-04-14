#ifndef _DBSCAN_CPU_H_
#define _DBSCAN_CPU_H_

#include <vector>
#include "datapoint.h"
#include "rtree.h"

#define NOISE -1

class DBSCAN {
    protected:
        int minPts; // Minimum number of points in a neighborhood to form a dense region
        double eps; // Maximum distance between two points to be considered neighbors

        int cluster_id; // Current cluster ID
        int *labels; // Array to store cluster labels for each point

        std::vector<DataPoint>* data; // pointer to the data points
        size_t data_size; // Size of the data array

        RTree<double, double, 3, double> tree; // RTree for spatial indexing
    private:
        // Function to calculate the distance between two points
        double getDist(DataPoint& p1, DataPoint& p2); 

        // Function to find all points within eps distance from the given point
        std::vector<int> regionQuery(int point, std::vector<DataPoint>& points); 

    public:
        DBSCAN(std::vector<DataPoint>& points, int minpts, double eps);
        ~DBSCAN();
        void run();
        void results();
};

#endif // _DBSCAN_CPU_H_