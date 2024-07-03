#include <iostream>
#include <vector>

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>


// __global__ void removeDuplicatesGPU(std::vector<unsigned long>& h_vec) {
//     int n = h_vec.size();

//     // Create a device vector and copy data from host to device
//     thrust::device_vector<unsigned long> d_vec = h_vec;

//     // Sort the vector
//     thrust::sort(d_vec.begin(), d_vec.end());

//     // Remove duplicates
//     thrust::device_vector<unsigned long>::iterator end = thrust::unique(d_vec.begin(), d_vec.end());

//     // Resize the vector to remove the undefined elements at the end
//     d_vec.erase(end, d_vec.end());

//     // Copy results back to host
//     h_vec.resize(d_vec.size());
//     thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());
// }

// __global__ void removeDuplicatesKernel(unsigned long* d_data, size_t num_elements) {
//     int tid = threadIdx.x + blockIdx.x * blockDim.x;
//     if (tid < num_elements - 1) {
//         if (d_data[tid] == d_data[tid + 1]) {
//             d_data[tid] = 0; // Marking duplicates
//         }
//     }
//     // Ensure all threads have completed marking duplicates before compacting
//     __syncthreads();

//     // Compact the data by shifting non-zero elements forward
//     if (tid < num_elements - 1) {
//         if (d_data[tid] == 0) {
//             for (int i = tid; i < num_elements - 1; ++i) {
//                 d_data[i] = d_data[i + 1];
//             }
//             d_data[num_elements - 1] = 0; // Set last element to zero
//         }
//     }
// }

void removeDuplicatesGPU(thrust::device_vector<unsigned long>& hash_vector) {
    // Sort the vector first
    thrust::sort(hash_vector.begin(), hash_vector.end());

    // Use thrust::unique to remove consecutive duplicates
    auto new_end = thrust::unique(hash_vector.begin(), hash_vector.end());

    // Resize the vector to remove the duplicates
    hash_vector.erase(new_end, hash_vector.end());
}