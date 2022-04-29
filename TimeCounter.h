#include <chrono>
decltype(std::chrono::system_clock::now()) startTime, endTime;
#define TIME_START startTime = std::chrono::system_clock::now();
#define TIME_END endTime = std::chrono::system_clock::now();\
    std::cout << "time:" << double(std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count())/1000000 << std::endl;
