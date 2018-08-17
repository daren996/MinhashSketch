#ifndef MINHASHSKETCH_UTILS_H
#define MINHASHSKETCH_UTILS_H

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <map>

//typedef unsigned int uint;
//typedef unsigned long ulong;
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

using namespace std;

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        std::cerr << "GPUassert: " << cudaGetErrorString(code) << ' ' << file << ':' << line << std::endl;
        if (abort)
            exit(code);
    }
}

namespace utils {
    int base2int(char base);
    void file_to_string(ifstream &file, string &info, string &sequence);
    uint pow_mod(long a, long b, long mod);
    bool isTextFile(const string& fileName);
}

std::ostream& bold_on(std::ostream& os);
std::ostream& bold_off(std::ostream& os);
std::ostream& uline_on(std::ostream& os);
std::ostream& uline_off(std::ostream& os);

#endif //MINHASHSKETCH_UTILS_H
