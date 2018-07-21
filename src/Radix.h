#ifndef MINHASHSKETCH_RADIX_H
#define MINHASHSKETCH_RADIX_H

#include <iostream>
#include <vector>
#include <cstdint>
#include <cmath>
#include <queue>
#include <thread>
#include <unordered_set>
#include "Hash.h"

namespace rs {

    typedef std::vector<uint64> List;

    signature genSig_single(int k, int m, string &sequence, vector<Hash> &hashes);
    signature genSig_multi(int k, int m, string &sequence, vector<Hash> &hashes);
    List rSort(int k, int m, string &sequence, uint64 begin, uint64 end, Hash &hash, List &list);
    void rSortLowestToHighest(List &list, uint64 begin, uint64 end);  // Sort starting from lowest digit.
    uint64 rMerge(List &list, uint64 begin1, uint64 end1, uint64 begin2, uint64 end2, int m);
}

#endif //MINHASHSKETCH_RADIX_H
