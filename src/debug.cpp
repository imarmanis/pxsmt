#include "debug.hpp"
#include <iostream>

bool debug_enabled = false;
bool trace_enabled = false;
bool prof_enabled = true;

void exit_msg(const char msg[])
{
    std::cerr << msg;
    exit(1);
}
