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

typedef unsigned long h_type;
typedef vector<vector<unsigned long>> signature;

class Hash {
private:
    long a;
    long b;
    long p;

    bool is_prime(long x);

    long generateNextPrime(long n);

public:
    Hash() = default;

    Hash(unsigned long u, int seed) {
        p = generateNextPrime(u);
        mt19937 rng(seed);
        uniform_int_distribution<h_type> distA(1, static_cast<h_type>(p - 1));
        uniform_int_distribution<h_type> distB(0, static_cast<h_type>(p - 1));
        a = distA(rng);
        b = distB(rng);
    }

    unsigned long operator()(unsigned long x) const;

    unsigned long operator()(uint64 *x, int k) const;
};

vector<Hash> generateHashes(int t, int seed);

void insertValue(unsigned long value, signature &s, unordered_set<unsigned long> &filter,
                 const vector<Hash> &hashes, vector<unsigned long> heap_max_v, int k);

signature generateSignature(int k, int m, const string &file, const vector<Hash> &hashes);

int computeSim(vector<unsigned long> v1, vector<unsigned long> v2);

double computeSim(const signature &sig1, const signature &sig2);

double computeSim(const string &file1, const string &file2, int k, int m, int t, int seed);

bool cmp(int a, int b);

#endif //MINHASHSKETCH_MINHASH_H
