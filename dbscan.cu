#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <vector>
#include <math.h>

// Structure to store 3D points
struct Point3D {
    float x, y, z;
};

// Constants for DBSCAN
const int NOISE = -1;
const int UNCLASSIFIED = -2;

// Kernel to calculate distance matrix between points
__global__ void calculateDistanceMatrix(Point3D* points, int numPoints, float* distMatrix, float eps) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (idx < numPoints && idy < numPoints) {
        float dx = points[idx].x - points[idy].x;
        float dy = points[idx].y - points[idy].y;
        float dz = points[idx].z - points[idy].z;
        
        float dist = sqrt(dx*dx + dy*dy + dz*dz);
        
        // Store 1 if points are within eps distance, 0 otherwise
        distMatrix[idx * numPoints + idy] = (dist <= eps) ? 1.0f : 0.0f;
    }
}

// Kernel to count neighbors for each point
__global__ void countNeighbors(float* distMatrix, int numPoints, int* neighborCounts) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < numPoints) {
        int count = 0;
        for (int j = 0; j < numPoints; j++) {
            if (distMatrix[idx * numPoints + j] > 0) {
                count++;
            }
        }
        neighborCounts[idx] = count;
    }
}

// Kernel to expand clusters
__global__ void expandCluster(float* distMatrix, int numPoints, int* labels, 
                            int* neighborCounts, int minPts, int currentCluster, 
                            bool* changed, int* borderPoints, int borderPointsCount) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < borderPointsCount) {
        int pointIdx = borderPoints[idx];
        
        // Skip already processed points
        if (labels[pointIdx] != currentCluster) {
            return;
        }
        
        for (int j = 0; j < numPoints; j++) {
            if (distMatrix[pointIdx * numPoints + j] > 0 && labels[j] == UNCLASSIFIED) {
                // Mark as part of the current cluster
                atomicExch(&labels[j], currentCluster);
                *changed = true;
                
                // If core point, add to border points for next iteration
                if (neighborCounts[j] >= minPts) {
                    // We'll handle this point in the next iteration
                }
            }
        }
    }
}

// Kernel to find initial core points and mark them with their cluster IDs
__global__ void findCorePoints(int* neighborCounts, int numPoints, int minPts, int* labels) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < numPoints) {
        if (neighborCounts[idx] >= minPts && labels[idx] == UNCLASSIFIED) {
            // This is a core point that hasn't been assigned to a cluster yet
        }
    }
}

// Host function to run DBSCAN
void dbscan(Point3D* h_points, int numPoints, float eps, int minPts, int* h_labels) {
    // Allocate device memory
    Point3D* d_points;
    float* d_distMatrix;
    int* d_labels;
    int* d_neighborCounts;
    
    cudaMalloc((void**)&d_points, numPoints * sizeof(Point3D));
    cudaMalloc((void**)&d_distMatrix, numPoints * numPoints * sizeof(float));
    cudaMalloc((void**)&d_labels, numPoints * sizeof(int));
    cudaMalloc((void**)&d_neighborCounts, numPoints * sizeof(int));
    
    // Initialize all points as UNCLASSIFIED
    for (int i = 0; i < numPoints; i++) {
        h_labels[i] = UNCLASSIFIED;
    }
    
    // Copy data to device
    cudaMemcpy(d_points, h_points, numPoints * sizeof(Point3D), cudaMemcpyHostToDevice);
    cudaMemcpy(d_labels, h_labels, numPoints * sizeof(int), cudaMemcpyHostToDevice);
    
    // Calculate distance matrix
    dim3 blockDim(16, 16);
    dim3 gridDim((numPoints + blockDim.x - 1) / blockDim.x, 
                (numPoints + blockDim.y - 1) / blockDim.y);
    
    calculateDistanceMatrix<<<gridDim, blockDim>>>(d_points, numPoints, d_distMatrix, eps);
    cudaDeviceSynchronize();
    
    // Count neighbors for each point
    int threadsPerBlock = 256;
    int blocksPerGrid = (numPoints + threadsPerBlock - 1) / threadsPerBlock;
    
    countNeighbors<<<blocksPerGrid, threadsPerBlock>>>(d_distMatrix, numPoints, d_neighborCounts);
    cudaDeviceSynchronize();
    
    // Copy neighbor counts back to host for cluster expansion
    int* h_neighborCounts = (int*)malloc(numPoints * sizeof(int));
    cudaMemcpy(h_neighborCounts, d_neighborCounts, numPoints * sizeof(int), cudaMemcpyDeviceToHost);
    
    // Main DBSCAN algorithm
    int currentCluster = 0;
    
    for (int i = 0; i < numPoints; i++) {
        if (h_labels[i] != UNCLASSIFIED) {
            continue;  // Skip already processed points
        }
        
        if (h_neighborCounts[i] < minPts) {
            h_labels[i] = NOISE;  // Mark as noise
            continue;
        }
        
        // New cluster found, expand it
        currentCluster++;
        h_labels[i] = currentCluster;
        
        // Use vectors to track border points that need processing
        std::vector<int> borderPoints;
        borderPoints.push_back(i);
        
        int borderIndex = 0;
        
        // Expand the cluster using both host and device
        while (borderIndex < borderPoints.size()) {
            int currentBorderSize = borderPoints.size() - borderIndex;
            int* d_borderPoints;
            bool* d_changed;
            bool h_changed = false;
            
            cudaMalloc((void**)&d_borderPoints, currentBorderSize * sizeof(int));
            cudaMalloc((void**)&d_changed, sizeof(bool));
            
            cudaMemcpy(d_borderPoints, &borderPoints[borderIndex], 
                    currentBorderSize * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(d_changed, &h_changed, sizeof(bool), cudaMemcpyHostToDevice);
            cudaMemcpy(d_labels, h_labels, numPoints * sizeof(int), cudaMemcpyHostToDevice);
            
            // Process current set of border points
            threadsPerBlock = 256;
            blocksPerGrid = (currentBorderSize + threadsPerBlock - 1) / threadsPerBlock;
            
            expandCluster<<<blocksPerGrid, threadsPerBlock>>>(
                d_distMatrix, numPoints, d_labels, d_neighborCounts, 
                minPts, currentCluster, d_changed, d_borderPoints, currentBorderSize);
            cudaDeviceSynchronize();
            
            // Get updated labels and changed flag
            cudaMemcpy(h_labels, d_labels, numPoints * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(&h_changed, d_changed, sizeof(bool), cudaMemcpyDeviceToHost);
            
            // Add new border points
            for (int j = 0; j < numPoints; j++) {
                if (h_labels[j] == currentCluster && 
                    std::find(borderPoints.begin(), borderPoints.end(), j) == borderPoints.end()) {
                    if (h_neighborCounts[j] >= minPts) {
                        borderPoints.push_back(j);
                    }
                }
            }
            
            borderIndex += currentBorderSize;
            
            cudaFree(d_borderPoints);
            cudaFree(d_changed);
            
            if (!h_changed) {
                break;  // No more expansion possible
            }
        }
    }
    
    // Clean up
    cudaFree(d_points);
    cudaFree(d_distMatrix);
    cudaFree(d_labels);
    cudaFree(d_neighborCounts);
    free(h_neighborCounts);
}

// Main function to run the DBSCAN clustering
int main(int argc, char** argv) {
    // Parse command line arguments (if any)
    float eps = 0.5f;  // Default epsilon value
    int minPts = 5;    // Default minimum points
    int numPoints = 10000; // Default number of points
    
    if (argc > 1) numPoints = atoi(argv[1]);
    if (argc > 2) eps = atof(argv[2]);
    if (argc > 3) minPts = atoi(argv[3]);
    
    printf("Running DBSCAN with params: numPoints=%d, eps=%.2f, minPts=%d\n", 
           numPoints, eps, minPts);
    
    // Allocate host memory for points and labels
    Point3D* h_points = (Point3D*)malloc(numPoints * sizeof(Point3D));
    int* h_labels = (int*)malloc(numPoints * sizeof(int));
    
    // Generate random points for testing
    srand(42);  // Fixed seed for reproducibility
    for (int i = 0; i < numPoints; i++) {
        h_points[i].x = ((float)rand() / RAND_MAX) * 100.0f;
        h_points[i].y = ((float)rand() / RAND_MAX) * 100.0f;
        h_points[i].z = ((float)rand() / RAND_MAX) * 100.0f;
    }
    
    // Record start time
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    
    // Run DBSCAN
    dbscan(h_points, numPoints, eps, minPts, h_labels);
    
    // Record end time
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    
    // Count clusters
    int maxCluster = 0;
    int noiseCount = 0;
    
    for (int i = 0; i < numPoints; i++) {
        if (h_labels[i] > maxCluster) {
            maxCluster = h_labels[i];
        }
        if (h_labels[i] == NOISE) {
            noiseCount++;
        }
    }
    
    printf("DBSCAN completed in %.2f ms\n", milliseconds);
    printf("Found %d clusters and %d noise points\n", maxCluster, noiseCount);
    
    // Clean up
    free(h_points);
    free(h_labels);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    
    return 0;
}