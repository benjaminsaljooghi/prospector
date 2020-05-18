#pragma once

#include <chrono>


inline std::chrono::high_resolution_clock::time_point time()
{

    return std::chrono::high_resolution_clock::now();
}

inline std::chrono::high_resolution_clock::time_point time(const std::chrono::high_resolution_clock::time_point& start, const char* message, const char* indent)
{
    auto curr = time();
    printf("%s%ldms %s\n", indent, std::chrono::duration_cast<std::chrono::milliseconds>(curr - start).count(), message);
    return curr;
}

inline std::chrono::high_resolution_clock::time_point time(const std::chrono::high_resolution_clock::time_point& start, const char* message)
{
	return time(start, message, "");
}


