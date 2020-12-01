#pragma once

#include <chrono>


inline std::chrono::high_resolution_clock::time_point time()
{

    return std::chrono::high_resolution_clock::now();
}

inline long long time_diff(const std::chrono::high_resolution_clock::time_point& start, const std::chrono::high_resolution_clock::time_point& curr)
{
    return std::chrono::duration_cast<std::chrono::milliseconds>(curr - start).count();
}

inline std::chrono::high_resolution_clock::time_point time(const std::chrono::high_resolution_clock::time_point& start, const char* message, const char* indent)
{
    auto curr = time();
    printf("%s%ldms %s\n", indent, time_diff(start, curr), message);
    return curr;
}

inline std::chrono::high_resolution_clock::time_point time(const std::chrono::high_resolution_clock::time_point& start, const char* message)
{
	return time(start, message, "");
}


