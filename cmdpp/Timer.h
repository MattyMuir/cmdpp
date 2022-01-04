#pragma once
#include <chrono>

class Timer
{
public:
    long long duration;
    Timer()
    {
        duration = 0;
        m_StartTimePoint = std::chrono::high_resolution_clock::now();
    }

    void Stop()
    {
        auto endTimePoint = std::chrono::high_resolution_clock::now();

        auto start = std::chrono::time_point_cast <std::chrono::microseconds> (m_StartTimePoint).time_since_epoch().count();
        auto end = std::chrono::time_point_cast <std::chrono::microseconds> (endTimePoint).time_since_epoch().count();

        duration = end - start;
    }
private:
    std::chrono::time_point < std::chrono::high_resolution_clock > m_StartTimePoint;
};