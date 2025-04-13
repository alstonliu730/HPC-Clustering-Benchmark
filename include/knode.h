#ifndef _KNODE_H_
#define _KNODE_H_

#include "datapoint.h"
#include <limits>

template <typename T>
class KNode {
    T *data; // Pointer to the data array
    int dim; // Dimension of the data point

    KNode<T> *left; // Pointer to the left child node
    KNode<T> *right; // Pointer to the right child node

    bool is_root; // Flag to indicate if this node is the root of the tree

    public:
        KNode(double *values, int dims, KNode<T> *l, KNode<T> *r, bool root)
            : data(new T[dim]), dim(dims), left(l), right(r), is_root(root) {}
        
        KNode(KNode<T> &&other)
            : data(other.data), dim(other.dim), 
            left(other.left), right(other.right),
            is_root(other.is_root) {
                other.data = nullptr; // Prevent double deletion
                other.left = nullptr;
                other.right = nullptr;
                other.dim = -1;
        }

        KNode():
            data(nullptr), dim(0), left(nullptr), right(nullptr), is_root(false) {}
        
        ~KNode() {
            if (is_root)
                delete[] data; // Free the allocated memory
            delete left;
            delete right;
        }

        /**
         * @brief the value of the data at the given index
         */
        T get_data(int index) const { return data[index]; }

        /**
         * @brief the number of dimensions of the data point
         */
        int get_dim() const { return dim; }

        /**
         * @brief Get pointer to the left child node
         */
        KNode<T>* get_left() const { return left; }

        /**
         * @brief Get pointer to the right child node
         */
        KNode<T>* get_right() const { return right; }

        /**
         * @brief check if this node is the root of the tree
         */
        bool is_root_node() const { return is_root; }
};  

#endif // KNODE_H