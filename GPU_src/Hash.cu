//
// Created by Darren on 2018/7/6.
//

#include "Hash.h"

bool Hash::is_prime(long long int x) {
    long long int i = 3;
    while (true) {
        long q = x / i;
        if (q < i) {
            return true;
        }
        if (x == q * i) {
            return false;
        }
        i += 2;
    }
}

// Tip: Save next prime of 2^62 to get it faster
long long int Hash::generateNextPrime(long long int n) {
    if (n <= 2) {
        return 2;
    }
    if (!(n & 1)) {
        ++n;
    }
    while (!is_prime(n)) {
        n += 2;
    }
    return n;
}

uint64 Hash::operator()(uint64 x) const {
    return (a * x + b) % p;
}

uint64 Hash::operator()(uint64 *x, int k) const {
    return SpookyHash::Hash64(x, k, b); // k is length of sequences in bytes
}

__device__
uint64 Hash::getValue_d(uint64 *x) {
    return (a * x[0] + b) % p;
}


vector<Hash> generateHashes(int t, int seed) {
    mt19937 rng(seed);
    uniform_int_distribution<int> distribution(0, INT_MAX);
    vector<Hash> hashes(t);
    for (int i = 0; i < t; ++i) {
        hashes[i] = Hash(LONG_MAX, distribution(rng));
    }
    return hashes;
}

int computeSim(vector<uint64> v1, vector<uint64> v2) {
    int i = 0, j = 0, count = 0;
    while (i < v1.size() && j < v2.size()) {
        if (v1[i] == v2[j]) {
            count += 1;
            i++;
            j++;
        } else if (v1[i] > v2[j])
            j++;
        else
            i++;
    }
    return count;
}

double computeSim(const signature &sig1, const signature &sig2) {
    int j = 0;
    for (int h = 0; h < sig1.size(); ++h) {
        j += computeSim(sig1[h], sig2[h]);
    }
    return double(j) / double(sig1.size() * sig1[0].size());
}


