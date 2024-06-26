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

__global__ void removeDuplicatesKernel(unsigned long* d_data, size_t num_elements) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid < num_elements - 1) {
        if (d_data[tid] == d_data[tid + 1]) {
            d_data[tid] = 0; // or any value to mark duplicates
        }
    }
}

void removeDuplicatesGPU(thrust::device_vector<unsigned long>& hash_vector) {
    size_t num_elements = hash_vector.size();
    unsigned long* d_data = thrust::raw_pointer_cast(hash_vector.data());

    // Sort the data to bring duplicates next to each other
    thrust::sort(hash_vector.begin(), hash_vector.end());

    // Remove duplicates using a kernel
    removeDuplicatesKernel<<<(num_elements + 255) / 256, 256>>>(d_data, num_elements);

    // Remove marked duplicates
    auto new_end = thrust::remove(hash_vector.begin(), hash_vector.end(), 0);
    hash_vector.erase(new_end, hash_vector.end());
}