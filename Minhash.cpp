//
// Created by Darren on 2018/7/6.
//

#include "Minhash.h"

bool Hash::is_prime(long x) {
    long i = 3;
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

long Hash::generateNextPrime(long n) {
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

h_type Hash::operator()(h_type k) const {
    return (h_type) (((a * k + b) % p) % m);
}


vector<Hash> generateHashes(int t, int seed) {
    mt19937 rng(seed);
    uniform_int_distribution<int> distribution(0, INT_MAX);
    vector<Hash> hashes(t);
    for (int i = 0; i < t; ++i) {
        hashes[i] = Hash(INT_MAX, INT_MAX, distribution(rng));
    }
    return hashes;
}

void insertValue(int value, signature &sig, unordered_set<int> &filter, const vector<Hash> &hashes) {
    if (filter.insert(value).second) {
        for (int h = 0; h < hashes.size(); ++h) {
            int p = hashes[h](static_cast<h_type>(value));
            sig[h].push_back(p);
            push_heap(sig[h].begin(), sig[h].end(), cmp);
            pop_heap(sig[h].begin(), sig[h].end(), cmp);
            sig[h].pop_back();
        }
    }
}

signature generateSignature(int k, int m, const string &sequence, const vector<Hash> &hashes) {
    signature sig(hashes.size(), vector<int>(m, INT_MAX));
    long s_index = 0; // pointer of current base
    for (int h = 0; h < hashes.size(); ++h) {
        make_heap(sig[h].begin(), sig[h].end(), cmp);
    }
    unordered_set<int> filter;
    hash<string> hashFunction; // transfer string to int by STL hash
    string shingle(static_cast<unsigned long>(k), ' ');
    for (int i = 0; i <= sequence.size() - k; ++i) {
        for (int j = 0; j < k; ++j) {
            shingle[j] = sequence[i + j];
        }
        insertValue(static_cast<int>(hashFunction(shingle)), sig, filter, hashes);
    }
    for (int h = 0; h < sig.size(); ++h) {
        sort_heap(sig[h].begin(), sig[h].end(), cmp);
    }
    return sig;
}

int computeSim(vector<int> v1, vector<int> v2) {
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

double computeSim(const string &sequence1, const string &sequence2, int k, int m, int t, int seed) {
    vector<Hash> hashes = generateHashes(t, seed);
    signature sig1 = generateSignature(k, m, sequence1, hashes);
    signature sig2 = generateSignature(k, m, sequence2, hashes);
    return computeSim(sig1, sig2);
}

bool cmp(int a, int b) {
    return a < b;
}

