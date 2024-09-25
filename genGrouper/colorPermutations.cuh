// File: color_permutations.cuh

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <cuda.h>
#include <iostream>
#include <vector>
#include <utility>

// Struct to represent a coloring with fixed-size array
struct Coloring {
    int* colors;
    int M;

    // Constructor to initialize colors array
    Coloring(int size) : M(size) {
        colors = new int[M];  // Dynamically allocate array of size M
    }

    // Destructor to clean up dynamically allocated memory
    ~Coloring() {
        delete[] colors;
    }

    // Copy constructor and assignment operator
    Coloring(const Coloring& other) : M(other.M) {
        colors = new int[M];
        std::copy(other.colors, other.colors + M, colors);
    }

    Coloring& operator=(const Coloring& other) {
        if (this != &other) {
            delete[] colors;
            M = other.M;
            colors = new int[M];
            std::copy(other.colors, other.colors + M, colors);
        }
        return *this;
    }

    // Comparison operator for sorting
    __host__ __device__
    bool operator<(const Coloring& other) const {
        for (int i = 0; i < M; ++i) {
            if (colors[i] < other.colors[i]) return true;
            if (colors[i] > other.colors[i]) return false;
        }
        return false;
    }

    // Equality operator for unique
    __host__ __device__
    bool operator==(const Coloring& other) const {
        for (int i = 0; i < M; ++i) {
            if (colors[i] != other.colors[i]) return false;
        }
        return true;
    }
};

// Function to apply automorphisms and collect unique colorings using GPU
std::vector<std::vector<int>> apply_edge_automorphisms_gpu(
    const std::vector<std::vector<int>>& colorings, 
    const std::vector<std::vector<std::pair<int, int>>>& automorphisms,
    const int M
);