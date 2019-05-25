#include <iostream>
#include <algorithm>
#include <string>
#include "../include/main.h"
#include "../include/util.h"

void time_stat(clock_t start_time, clock_t end_time) {
    /* how to time a parallel program? */
    clock_t duration = end_time - start_time;
    size_t duration_sec = duration / CLOCKS_PER_SEC;
    size_t duration_min = duration_sec / 60;
    size_t duration_hr  = duration_min / 60;
    duration_sec = duration_sec % 60;
    duration_min = duration_min % 60;
    std::cout << "EZNMA> Elapsed time: "
              << duration_hr  << " HR "
              << duration_min << " MIN "
              << duration_sec << " SEC." << std::endl;
}

std::string real2str(real inp_num) {
    auto tmp_str = std::to_string(inp_num);
    auto dot_pos = tmp_str.find(".");
    auto result = tmp_str.substr(0, dot_pos+2);
    return result;
}