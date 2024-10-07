#ifndef COLORING_CUH
#define COLORING_CUH

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <iostream>
#include <vector>

// Constants
#define M 10
#define PADDING -1

// Struct to represent a coloring with fixed-size array
struct Coloring {
    int colors[M];

    // Comparison operator for sorting
    __host__ __device__
    bool operator<(const Coloring& other) const {
        for (int i = 0; i < M; ++i) {
            if (colors[i] < other.colors[i]) return true;
            if (colors[i] > other.colors[i]) return false;
        }
        return false;
    }

    // Equality operator for uniqueness check
    __host__ __device__
    bool operator==(const Coloring& other) const {
        for (int i = 0; i < M; ++i) {
            if (colors[i] != other.colors[i]) return false;
        }
        return true;
    }
};

// CUDA kernel function declaration
__global__ void find_canonical_colorings_kernel(
    const int* colorings,              // N x M
    const int* automorphisms,          // K x M (source indices)
    Coloring* canonical_colorings,     // N x M
    int N,                             // Number of colorings
    int K,                             // Number of automorphisms
    int M_edges                        // Number of edges
);

// Function declaration for applying automorphisms and collecting unique colorings using GPU
std::vector<std::vector<int>> apply_edge_automorphisms_gpu(
    const std::vector<std::vector<int>>& colorings,
    const std::vector<std::vector<std::pair<int, int>>>& automorphisms
);

#endif // COLORING_CUH
