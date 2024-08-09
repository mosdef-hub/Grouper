#include <iostream>
#include <unordered_map>
#include <vector>
#include <sys/resource.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>




void logMemoryUsage(const std::string& message) {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    std::cout << message << " - Memory usage: " << usage.ru_maxrss << " kilobytes" << std::endl;
}

template<typename T>
size_t getHostVectorMemoryUsage(const thrust::host_vector<T>& vec) {
    return vec.capacity() * sizeof(T);
}

template<typename T>
size_t getVectorMemoryUsage(const std::vector<T>& vec) {
    return vec.capacity() * sizeof(T);
}

// Helper function to estimate memory usage of unordered maps
template<typename K, typename V>
size_t getUnorderedMapMemoryUsage(const std::unordered_map<K, V>& umap) {
    return umap.size() * (sizeof(K) + sizeof(V)) + sizeof(umap);
}

// Helper function to estimate memory usage of unordered sets
template<typename T>
size_t getUnorderedSetMemoryUsage(const std::unordered_set<T>& uset) {
    return uset.size() * sizeof(T) + sizeof(uset);
}

class Logger {
public:
    template <typename T>
    void logSize(const std::string& name, const T& var) {
        sizes[name] = sizeof(var);
        types[name] = typeid(var).name();
    }

    void printSizes() const {
        for (const auto& entry : sizes) {
            std::cout << "Size of " << entry.first << " (" << types.at(entry.first) << "): " << entry.second << " bytes" << std::endl;
        }
    }

private:
    std::unordered_map<std::string, size_t> sizes;
    std::unordered_map<std::string, std::string> types;
};