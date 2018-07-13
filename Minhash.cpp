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

// Tip: Save next prime of 2^62 to get it faster
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

unsigned long Hash::operator()(unsigned long x) const {
    return (a * x + b) % p;
}
unsigned long Hash::operator()(uint64 *x, int k) const {
    return SpookyHash::Hash64(x, k, b);
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

void insertValue(uint64 *value, signature &sig, unordered_set<uint64> &filter,
                 const vector<Hash> &hashes, vector<uint64> heap_max_v, int k) {
    if (filter.insert(value[0]).second) {
        for (int h = 0; h < hashes.size(); ++h) {
            uint64 p = hashes[h](value, k*2);
            if (p < heap_max_v[h]) {
                sig[h].push_back(p);
                push_heap(sig[h].begin(), sig[h].end(), cmp);
                pop_heap(sig[h].begin(), sig[h].end(), cmp);
                sig[h].pop_back();
            }
        }
    }
}

// Haven't written what to do if k is greater than 32
signature generateSignature(int k, int m, const string &sequence, const vector<Hash> &hashes) {
    uint64 s_index = 0; // pointer of current base
    uint64 cur_seq[k / 32 + 1]; // current sub-sequence
    for (int i = 0; i < k / 32 + 1; ++i) {
        cur_seq[i] = 0;
    }
    unordered_set<uint64> filter;
    signature sig(hashes.size(), vector<uint64>(m, UINT64_MAX ));
    vector<uint64> heap_max_v(m, UINT64_MAX );
    for (int h = 0; h < hashes.size(); ++h) {
        make_heap(sig[h].begin(), sig[h].end(), cmp);
    }
    if (k < 32) {
        for (; s_index < k; ++s_index) {
            cur_seq[0] = (cur_seq[0] << 2) % (uint64)pow(4, k) + utils::base2int(sequence[s_index]);
        }
        insertValue(cur_seq, sig, filter, hashes, heap_max_v, k);
        for (; s_index < sequence.size(); ++s_index) {
            cur_seq[0] = (cur_seq[0] << 2) % (uint64)pow(4, k) + utils::base2int(sequence[s_index]);
            insertValue(cur_seq, sig, filter, hashes, heap_max_v, k);
        }
    } else {
        for (; s_index < k; ++s_index) {
            cur_seq[s_index / 32] = (cur_seq[s_index / 32] << 2) % UINT64_MAX + utils::base2int(sequence[s_index]);
        }
        insertValue(cur_seq, sig, filter, hashes, heap_max_v, k);

        for (; s_index < sequence.size(); ++s_index) {
            for (int i = 0; i < k / 32 - 1; ++i) {
                cur_seq[i] = (cur_seq[i] << 2) + (cur_seq[i + 1] >> 62);
            }
            cur_seq[k / 32 - 1] = (cur_seq[k / 32 - 1] << 2) + (cur_seq[k / 32] >> ((k % 32) * 2));
            cur_seq[k / 32] = (cur_seq[k / 32] << 2) % (uint64)pow(4, k % 32) + utils::base2int(sequence[s_index]);
            insertValue(cur_seq, sig, filter, hashes, heap_max_v, k);
        }
    }
    for (int h = 0; h < sig.size(); ++h) {
        sort_heap(sig[h].begin(), sig[h].end(), cmp);
    }
    return sig;
}

int computeSim(vector<unsigned long> v1, vector<unsigned long> v2) {
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

