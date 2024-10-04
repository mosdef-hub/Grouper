// File: color_permutations.cu

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <cuda.h>
#include <iostream>
#include <vector>
#include <utility>

// Define the number of edges (fixed for simplicity)
#define M 10
#define PADDING -1

// Struct to represent a coloring with fixed-size array
struct Coloring {
    int colors[M];

    // Comparison operator for sorting
    __host__ __device__
    bool operator<(const Coloring& other) const {
        for(int i = 0; i < M; ++i){
            if(colors[i] < other.colors[i]) return true;
            if(colors[i] > other.colors[i]) return false;
        }
        return false;
    }

    // Equality operator for unique
    __host__ __device__
    bool operator==(const Coloring& other) const {
        for(int i = 0; i < M; ++i){
            if(colors[i] != other.colors[i]) return false;
        }
        return true;
    }
};

// CUDA kernel to find canonical colorings
__global__ void find_canonical_colorings_kernel(
    const int* colorings,              // N x M
    const int* automorphisms,          // K x M (source indices)
    Coloring* canonical_colorings,     // N x M
    int N,                              // Number of colorings
    int K,                              // Number of automorphisms
    int M_edges                         // Number of edges
){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= N) return;

    // Initialize min_coloring with the original coloring
    Coloring min_coloring;
    for(int i = 0; i < M_edges; ++i){
        min_coloring.colors[i] = colorings[idx * M_edges + i];
    }

    // Iterate over all automorphisms to find the minimal (canonical) coloring
    for(int k = 0; k < K; ++k){
        // Apply automorphism k
        Coloring permuted_coloring;
        for(int i = 0; i < M_edges; ++i){
            int source_idx = automorphisms[k * M_edges + i];
            permuted_coloring.colors[i] = colorings[idx * M_edges + source_idx];
        }

        // Compare permuted_coloring with min_coloring
        bool is_less = false;
        bool is_greater = false;
        for(int i = 0; i < M_edges; ++i){
            if(permuted_coloring.colors[i] < min_coloring.colors[i]){
                is_less = true;
                break;
            }
            if(permuted_coloring.colors[i] > min_coloring.colors[i]){
                is_greater = true;
                break;
            }
        }

        if(is_less){
            min_coloring = permuted_coloring;
        }
    }

    // Assign the minimal (canonical) coloring
    canonical_colorings[idx] = min_coloring;
}

// Function to apply automorphisms and collect unique colorings using GPU
std::vector<std::vector<int>> apply_edge_automorphisms_gpu(
    const std::vector<std::vector<int>>& colorings, 
    const std::vector<std::vector<std::pair<int, int>>>& automorphisms
){
    size_t N = colorings.size(); // Number of colorings
    size_t K = automorphisms.size(); // Number of automorphisms

    // Flatten colorings: N x M
    thrust::host_vector<int> h_colorings_flat;
    h_colorings_flat.reserve(N * M);
    for(const auto& coloring : colorings){
        for(int i = 0; i < M; ++i){
            int color = i < coloring.size() ? coloring[i] : PADDING;
            h_colorings_flat.push_back(color);
        }
    }
    thrust::device_vector<int> d_colorings_flat = h_colorings_flat;

    // Flatten automorphisms: K x M (store source indices)
    thrust::host_vector<int> h_automorphisms_flat;
    h_automorphisms_flat.reserve(K * M);
    for(const auto& automorphism : automorphisms){
        for(int i = 0; i < M; ++i){
            if (i < automorphism.size()){
                h_automorphisms_flat.push_back(automorphism[i].first);
            } else {
                h_automorphisms_flat.push_back(PADDING);
            }
        }
    }
    thrust::device_vector<int> d_automorphisms_flat = h_automorphisms_flat;

    // Allocate space for canonical colorings: N x M
    thrust::device_vector<Coloring> d_canonical_colorings(N);

    // Launch canonical mapping kernel
    int threads_per_block = 256;
    int blocks = (N + threads_per_block - 1) / threads_per_block;
    find_canonical_colorings_kernel<<<blocks, threads_per_block>>>(
        thrust::raw_pointer_cast(d_colorings_flat.data()),
        thrust::raw_pointer_cast(d_automorphisms_flat.data()),
        thrust::raw_pointer_cast(d_canonical_colorings.data()),
        N, K, M
    );

    // Error checking
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess){
        std::cerr << "CUDA kernel launch error: " << cudaGetErrorString(err) << std::endl;
        exit(EXIT_FAILURE);
    }
    cudaDeviceSynchronize();

    // Sort the canonical colorings lexicographically
    thrust::sort(d_canonical_colorings.begin(), d_canonical_colorings.end());

    // Remove duplicates
    auto new_end = thrust::unique(d_canonical_colorings.begin(), d_canonical_colorings.end());

    // Erase the duplicates
    d_canonical_colorings.erase(new_end, d_canonical_colorings.end());

    // Transfer unique colorings back to host
    thrust::host_vector<Coloring> h_unique_colorings_struct = d_canonical_colorings;

    // Convert structs back to vector of vectors
    std::vector<std::vector<int>> unique_colorings;
    unique_colorings.reserve(h_unique_colorings_struct.size());
    for(const auto& coloring_struct : h_unique_colorings_struct){
        std::vector<int> coloring;
        coloring.reserve(M);
        for(int i = 0; i < M; ++i){
            coloring.push_back(coloring_struct.colors[i]);
        }
        unique_colorings.push_back(coloring);
    }

    return unique_colorings;
}

// int main(){

//     // Define colorings for a triangle graph with 2 colors (0 and 1)
//     std::vector<std::vector<int>> tri_colorings = {
//         {0, 0, 0},
//         {0, 0, 1},
//         {0, 1, 0},
//         {0, 1, 1},
//         {1, 0, 0},
//         {1, 0, 1},
//         {1, 1, 0},
//         {1, 1, 1}
//     };

//     // Define automorphisms (identity, rotations, and reflections for triangle)
//     std::vector<std::vector<std::pair<int, int>>> tri_automorphisms = {
//         { {0,0}, {1,1}, {2,2} }, // Identity
//         { {1,1}, {2,2}, {0,0} }, // Rotation 1
//         { {2,2}, {0,0}, {1,1} }, // Rotation 2
//         { {0,0}, {2,2}, {1,1} }, // Reflection 1
//         { {1,1}, {0,0}, {2,2} }, // Reflection 2
//         { {2,2}, {1,1}, {0,0} }  // Reflection 3
//     };

//     // Generate all possible colorings for 6 edges with 2 colors
//     std::vector<std::vector<int>> hex_colorings;
//     for(int i = 0; i < 64; ++i){
//         std::vector<int> coloring(M);
//         for(int j = 0; j < M; ++j){
//             coloring[j] = (i >> j) & 1; // Extract the j-th bit
//         }
//         coloring[M -1] = PADDING;
//         hex_colorings.push_back(coloring);
//     }
//     // Define automorphisms (identity, rotations, and reflections for hexagon)
//     std::vector<std::vector<std::pair<int, int>>> hex_automorphisms = {
//         // Identity
//         { {0,0}, {1,1}, {2,2}, {3,3}, {4,4}, {5,5} },
//         // 60° Rotation
//         { {5,0}, {0,1}, {1,2}, {2,3}, {3,4}, {4,5} },
//         // 120° Rotation
//         { {4,0}, {5,1}, {0,2}, {1,3}, {2,4}, {3,5} },
//         // 180° Rotation
//         { {3,0}, {4,1}, {5,2}, {0,3}, {1,4}, {2,5} },
//         // 240° Rotation
//         { {2,0}, {3,1}, {4,2}, {5,3}, {0,4}, {1,5} },
//         // 300° Rotation
//         { {1,0}, {2,1}, {3,2}, {4,3}, {5,4}, {0,5} },
//         // Reflection over Axis through Edge 0-3
//         { {0,0}, {1,5}, {2,4}, {3,3}, {4,2}, {5,1} },
//         // Reflection over Axis through Edge 1-4
//         { {5,1}, {1,1}, {0,5}, {3,3}, {2,4}, {4,2} },
//         // Reflection over Axis through Edge 2-5
//         { {4,2}, {3,4}, {2,2}, {1,5}, {0,0}, {5,1} },
//         // Reflection over Axis through Vertex 0 and Vertex 3
//         { {0,0}, {5,1}, {4,2}, {3,3}, {2,4}, {1,5} },
//         // Reflection over Axis through Vertex 1 and Vertex 4
//         { {1,1}, {0,0}, {5,2}, {4,3}, {3,4}, {2,5} },
//         // Reflection over Axis through Vertex 2 and Vertex 5
//         { {2,2}, {1,1}, {0,0}, {5,5}, {4,4}, {3,3} }
//     };

//     // Apply edge automorphisms using GPU
//     std::vector<std::vector<int>> tri_unique_colorings = apply_edge_automorphisms_gpu(tri_colorings, tri_automorphisms);

//     // Print unique colorings
//     std::cout << "Unique Colorings:\n";
//     for(const auto& coloring : tri_unique_colorings){
//         for(const auto& color : coloring){
//             std::cout << color << " ";
//         }
//         std::cout << "\n";
//     }
//     std::cout << "Number of unique colorings: " << tri_unique_colorings.size() << "\n";

//     // Apply edge automorphisms using GPU
//     std::vector<std::vector<int>> hex_unique_colorings = apply_edge_automorphisms_gpu(hex_colorings, hex_automorphisms);

//         // Print unique colorings
//     std::cout << "Unique Colorings:\n";
//     for(const auto& coloring : hex_unique_colorings){
//         for(const auto& color : coloring){
//             std::cout << color << " ";
//         }
//         std::cout << "\n";
//     }
//     std::cout << std::flush;

//     // Print the number of unique colorings
//     std::cout<< "Number of unique colorings: " << hex_unique_colorings.size() << "\n";

//     return 0;
// }
