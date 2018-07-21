//
// Created by Darren on 2018/7/6.
//

#ifndef MINHASHSKETCH_MINHASH_H
#define MINHASHSKETCH_MINHASH_H

#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <random>
#include <climits>
#include <unordered_set>
#include "Utils.h"
#include "SpookyV2.h"

using namespace std;

typedef vector<vector<uint64>> signature;

class Hash {
private:
    uint64 a;
    uint64 b;
    uint64 p;

    bool is_prime(long long int x);

    long long int generateNextPrime(long long int n);

public:
    Hash() = default;

    Hash(uint64 u, int seed) {
//        p = generateNextPrime(static_cast<long long int>(u));
        p = (uint64)13835058055282163729;
        mt19937 rng(seed);
        uniform_int_distribution<uint64> distA(1, p - 1);
        uniform_int_distribution<uint64> distB(0, p - 1);
        a = distA(rng);
        b = distB(rng);
    }

    uint64 operator()(uint64 x) const;

    uint64 operator()(uint64 *x, int k) const;
};

vector<Hash> generateHashes(int t, int seed);

void insertValue(uint64 *value, signature &s, unordered_set<uint64> *filter,
                 const vector<Hash> &hashes, vector<uint64> heap_max_v, int k);

signature generateSignature(int k, int m, const string &file, const vector<Hash> &hashes);

int computeSim(vector<uint64> v1, vector<uint64> v2);

double computeSim(const signature &sig1, const signature &sig2);

double computeSim(const string &file1, const string &file2, int k, int m, int t, int seed);

bool cmp(int a, int b);

#endif //MINHASHSKETCH_MINHASH_H
