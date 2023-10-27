#pragma once

#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

std::string timestamp() {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);

    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X");
    return ss.str();
}

class logger {
  public:
    static std::ostream &stream() { return std::cout; }
};

std::ostream &log() { return (logger::stream() << "[" << timestamp()) << "] "; }
std::ostream &log(std::string pre) { return (logger::stream() << "[" << timestamp()) << "] [" << pre << "] "; }

class timer {

  public:
    timer() { start = std::chrono::system_clock::now(); }

    std::string duration() {
        auto total = std::chrono::system_clock::now() - start;
        const auto mins = std::chrono::duration_cast<std::chrono::minutes>(total);
        const auto secs = std::chrono::duration_cast<std::chrono::seconds>(total - mins);
        std::stringstream ss;
        ss << mins.count() << "m " << secs.count() << "s";
        return ss.str();
    }

  private:
    std::chrono::time_point<std::chrono::system_clock> start;
};
