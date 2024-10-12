#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <unordered_map>
#include <vector>
#include <sys/resource.h>
// #include <thrust/host_vector.h>
// #include <thrust/device_vector.h>
#include <fstream>
#include <string>
#include <ctime>
#include <mutex>




// void logMemoryUsage(const std::string& message) {
//     struct rusage usage;
//     getrusage(RUSAGE_SELF, &usage);
//     std::cout << message << " - Memory usage: " << usage.ru_maxrss << " kilobytes" << std::endl;
// }

// template<typename T>
// size_t getHostVectorMemoryUsage(const thrust::host_vector<T>& vec) {
//     return vec.capacity() * sizeof(T);
// }

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
// template<typename T>
// size_t getUnorderedSetMemoryUsage(const std::unordered_set<T>& uset) {
//     return uset.size() * sizeof(T) + sizeof(uset);
// }

class Logger {
public:
    enum LogLevel { DEBUG, INFO, WARNING, ERROR };

    // Singleton pattern for global access to the logger instance
    static Logger& getInstance() {
        static Logger instance;
        return instance;
    }

    // Set log level for filtering messages
    void setLogLevel(LogLevel level) {
        logLevel = level;
    }

    // Log to file and console
    void log(const std::string& message, LogLevel level) {
        std::lock_guard<std::mutex> guard(logMutex);
        if (level >= logLevel) {
            std::string levelStr = logLevelToString(level);
            std::string timeStr = currentDateTime();

            // Write to console
            std::cout << "[" << timeStr << "] " << levelStr << ": " << message << std::endl;

            // Write to file if file logging is enabled
            if (logFile.is_open()) {
                logFile << "[" << timeStr << "] " << levelStr << ": " << message << std::endl;
                logFile.flush();  // Ensure the log gets written to the file immediately
            }
        }
    }


    // Enable logging to file
    bool enableFileLogging(const std::string& filename) {
        std::lock_guard<std::mutex> guard(logMutex);
        logFile.open(filename, std::ios::out | std::ios::app);
        return logFile.is_open();
    }

    // Disable file logging
    void disableFileLogging() {
        std::lock_guard<std::mutex> guard(logMutex);
        if (logFile.is_open()) {
            logFile.close();
        }
    }

private:
    LogLevel logLevel = INFO; // Default log level
    std::ofstream logFile;
    std::mutex logMutex;

    Logger() = default;
    ~Logger() {
        if (logFile.is_open()) {
            logFile.close();
        }
    }

    // Delete copy constructor and assignment operator to enforce singleton
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

    // Convert log level to string
    std::string logLevelToString(LogLevel level) {
        switch (level) {
            case DEBUG: return "DEBUG";
            case INFO: return "INFO";
            case WARNING: return "WARNING";
            case ERROR: return "ERROR";
            default: return "UNKNOWN";
        }
    }

    // Get current date and time as a string
    std::string currentDateTime() {
        time_t now = time(0);
        tm* ltm = localtime(&now);
        char buffer[80];
        strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", ltm);
        return std::string(buffer);
    }
};

extern Logger& logger;

#endif // LOGGER_H
