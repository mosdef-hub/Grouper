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

// CUDA kernel remains the same
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
            if(source_idx == PADDING){
                permuted_coloring.colors[i] = PADDING;
            } else {
                permuted_coloring.colors[i] = colorings[idx * M_edges + source_idx];
            }
        }

        // Compare permuted_coloring with min_coloring
        bool is_less = false;
        for(int i = 0; i < M_edges; ++i){
            if(permuted_coloring.colors[i] < min_coloring.colors[i]){
                is_less = true;
                break;
            }
            if(permuted_coloring.colors[i] > min_coloring.colors[i]){
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
    size_t total_N = colorings.size(); // Total number of colorings
    size_t K = automorphisms.size(); // Number of automorphisms

    // std::cout << "Total number of colorings (N): " << total_N << std::endl;
    // std::cout << "Number of automorphisms (K): " << K << std::endl;

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

    // Determine available GPU memory
    size_t free_mem, total_mem;
    cudaError_t cudaStatus = cudaMemGetInfo(&free_mem, &total_mem);
    if(cudaStatus != cudaSuccess){
        std::cerr << "Failed to get GPU memory info: " << cudaGetErrorString(cudaStatus) << std::endl;
        exit(EXIT_FAILURE);
    }
    // std::cout << "Available GPU memory: " << free_mem / (1024 * 1024) << " MB out of " << total_mem / (1024 * 1024) << " MB\n";

    // Calculate memory required for each coloring batch
    size_t bytes_per_coloring = M * sizeof(int); // Original colorings
    size_t bytes_per_canonical = sizeof(Coloring); // Canonical colorings
    size_t bytes_automorphisms = K * M * sizeof(int);

    // Estimate maximum N that fits into free memory
    // Additionally, account for other allocations and some buffer (e.g., 10%)
    size_t buffer = free_mem / 10;
    size_t usable_mem = free_mem - buffer;

    // Each batch requires:
    // - colorings_flat: N_batch * M * sizeof(int)
    // - canonical_colorings: N_batch * sizeof(Coloring)
    // Total per batch: N_batch * (M * sizeof(int) + sizeof(Coloring))
    size_t memory_per_coloring = M * sizeof(int) + sizeof(Coloring);
    size_t max_N_per_batch = usable_mem / memory_per_coloring;

    if(max_N_per_batch == 0){
        std::cerr << "Not enough GPU memory to process even a single coloring.\n";
        exit(EXIT_FAILURE);
    }

    // std::cout << "Processing colorings in batches of up to " << max_N_per_batch << " colorings.\n";

    std::vector<std::vector<int>> unique_colorings_all;

    size_t processed_N = 0;
    while(processed_N < total_N){
        size_t current_batch_size = std::min(max_N_per_batch, total_N - processed_N);
        // std::cout << "Processing batch: " << processed_N << " to " << (processed_N + current_batch_size -1) << "\n";

        // Prepare the current batch of colorings
        thrust::host_vector<int> h_colorings_flat;
        h_colorings_flat.reserve(current_batch_size * M);
        for(size_t i = processed_N; i < processed_N + current_batch_size; ++i){
            for(int j = 0; j < M; ++j){
                int color = j < colorings[i].size() ? colorings[i][j] : PADDING;
                h_colorings_flat.push_back(color);
            }
        }
        thrust::device_vector<int> d_colorings_flat = h_colorings_flat;

        // Allocate space for canonical colorings: current_batch_size x M
        thrust::device_vector<Coloring> d_canonical_colorings(current_batch_size);

        // Launch canonical mapping kernel
        int threads_per_block = 256;
        int blocks = (current_batch_size + threads_per_block - 1) / threads_per_block;

        // std::cout << "Launching kernel with " << blocks << " blocks of " << threads_per_block << " threads.\n";

        find_canonical_colorings_kernel<<<blocks, threads_per_block>>>(
            thrust::raw_pointer_cast(d_colorings_flat.data()),
            thrust::raw_pointer_cast(d_automorphisms_flat.data()),
            thrust::raw_pointer_cast(d_canonical_colorings.data()),
            current_batch_size, K, M
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

        // Remove duplicates within the batch
        auto new_end = thrust::unique(d_canonical_colorings.begin(), d_canonical_colorings.end());

        // Erase the duplicates
        d_canonical_colorings.erase(new_end, d_canonical_colorings.end());

        // Transfer unique colorings back to host
        thrust::host_vector<Coloring> h_unique_colorings_struct = d_canonical_colorings;

        // Convert structs back to vector of vectors and add to the global list
        for(const auto& coloring_struct : h_unique_colorings_struct){
            std::vector<int> coloring;
            coloring.reserve(M);
            for(int i = 0; i < M; ++i){
                coloring.push_back(coloring_struct.colors[i]);
            }
            unique_colorings_all.push_back(coloring);
        }

        processed_N += current_batch_size;
    }

    // std::cout << "Total unique colorings before global deduplication: " << unique_colorings_all.size() << "\n";

    // At this point, unique_colorings_all contains unique colorings per batch.
    // Now, transfer all unique colorings to the GPU again to perform a global unique operation.

    // Flatten unique_colorings_all
    thrust::host_vector<int> h_all_unique_flat;
    h_all_unique_flat.reserve(unique_colorings_all.size() * M);
    for(const auto& coloring : unique_colorings_all){
        for(int i = 0; i < M; ++i){
            h_all_unique_flat.push_back(coloring[i]);
        }
    }

    thrust::device_vector<int> d_all_unique_flat = h_all_unique_flat;

    // Convert to Coloring structs
    struct ColoringComparator {
        __host__ __device__
        bool operator()(const Coloring& a, const Coloring& b) const {
            return a < b;
        }
    };

    size_t total_unique = unique_colorings_all.size();
    thrust::device_vector<Coloring> d_all_unique_colorings(total_unique);

    // Populate d_all_unique_colorings
    // This step can be optimized, but for simplicity, we'll do it on the host and transfer
    thrust::host_vector<Coloring> h_all_unique_colorings_host;
    h_all_unique_colorings_host.reserve(total_unique);
    for(const auto& coloring : unique_colorings_all){
        Coloring c;
        for(int i = 0; i < M; ++i){
            c.colors[i] = coloring[i];
        }
        h_all_unique_colorings_host.push_back(c);
    }
    d_all_unique_colorings = h_all_unique_colorings_host;

    // Sort all unique colorings globally
    thrust::sort(d_all_unique_colorings.begin(), d_all_unique_colorings.end());

    // Remove duplicates globally
    auto global_new_end = thrust::unique(d_all_unique_colorings.begin(), d_all_unique_colorings.end());

    // Erase the duplicates
    d_all_unique_colorings.erase(global_new_end, d_all_unique_colorings.end());

    // Transfer unique colorings back to host
    thrust::host_vector<Coloring> h_final_unique_colorings_struct = d_all_unique_colorings;

    // Convert structs back to vector of vectors
    std::vector<std::vector<int>> final_unique_colorings;
    final_unique_colorings.reserve(h_final_unique_colorings_struct.size());
    for(const auto& coloring_struct : h_final_unique_colorings_struct){
        std::vector<int> coloring;
        coloring.reserve(M);
        for(int i = 0; i < M; ++i){
            coloring.push_back(coloring_struct.colors[i]);
        }
        final_unique_colorings.push_back(coloring);
    }

    // std::cout << "Total unique colorings after global deduplication: " << final_unique_colorings.size() << "\n";

    return final_unique_colorings;
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