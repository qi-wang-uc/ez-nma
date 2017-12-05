#include <iostream>
#include <algorithm>
#include <string>
#include "../include/main.h"
#include "../include/util.h"

std::string real2str(real inp_num) {
    auto tmp_str = std::to_string(inp_num);
    auto dot_pos = tmp_str.find(".");
    auto result = tmp_str.substr(0, dot_pos+2);
    return result;
}