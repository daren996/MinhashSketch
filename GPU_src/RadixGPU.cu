#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <cstdint>
#include <cmath>
#include <queue>
#include <thread>
#include <unordered_set>
#include <cuda_runtime.h>
#include "Hash.h"
#include "SpookyV2_d.h"
#include <cub/cub.cuh>

using namespace std;

//#define BLOCKS_NUM 5
//#define BLOCK_THREADS 128
//#define ITEMS_PER_THREAD 4
//#define CEILING_DIVIDE(X, Y) (1 + (((X) - 1) / (Y)))
#define CHECK(res) if(res!=cudaSuccess){printf("CHECK ERROR!\n");exit(-1);}

__device__
int base2int(char base) {
    switch (base) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            return -1;
    }
}

__device__ uint64 getHashValue (uint64 *x, uint64 b, int k) {
    return SpookyHash_d::Hash64(x, (k / 32 + 1) * 8, b);
}

/*
 * Get hash value list on blocks
*/
template<int BLOCK_THREADS, int ITEMS_PER_THREAD>
__device__ void BlockGetHashValues (const int k, char *dna_d, int thread_offset, uint64 *thread_dataR, int hash_b) {
    int list_index = 0;
    int dna_index = thread_offset;
    bool isEnd = 0;
    uint64 *cur_seq = new uint64[k / 32 + 1]; // current sub-sequence
    for (int i = 0; i < k / 32 + 1; ++i)
        cur_seq[i] = 0;
    if (k < 32) {
        for (; dna_index < thread_offset + k - 1; ++dna_index) {
            if (dna_d[dna_index] == 'S')
                isEnd = 1;
            if (base2int(dna_d[dna_index]) != -1)
                cur_seq[0] = (cur_seq[0] << 2) % ((uint64) 1 << (2 * k)) + base2int(dna_d[dna_index]);
        }
        for (; dna_index < thread_offset + ITEMS_PER_THREAD + k - 1; ++dna_index) {
            if (dna_d[dna_index] == 'S')
                isEnd = 1;
            if (isEnd)
                thread_dataR[list_index++] = UINT64_MAX;
            else {
                if (base2int(dna_d[dna_index]) != -1)
                    cur_seq[0] = (cur_seq[0] << 2) % ((uint64) 1 << (2 * k)) + base2int(dna_d[dna_index]);
                thread_dataR[list_index++] = getHashValue(cur_seq, hash_b, k);
            }
        }
    } else {
        for (; dna_index < thread_offset + k; ++dna_index) {
            if (dna_d[dna_index] == 'S')
                isEnd = 1;
            if (base2int(dna_d[dna_index]) != -1)
                cur_seq[dna_index / 32] =
                        (cur_seq[dna_index / 32] << 2) % UINT64_MAX + base2int(dna_d[dna_index]);
        }
        if (isEnd)
            thread_dataR[list_index++] = UINT64_MAX;
        else
            thread_dataR[list_index++] = getHashValue(cur_seq, hash_b, k);
        for (; dna_index < thread_offset + ITEMS_PER_THREAD + k - 1; ++dna_index) {
            for (int j = 0; j < k / 32 - 1; ++j) {
                cur_seq[j] = (cur_seq[j] << 2) + (cur_seq[j + 1] >> 62);
            }
            cur_seq[k / 32 - 1] = (cur_seq[k / 32 - 1] << 2) + (cur_seq[k / 32] >> ((k % 32) * 2 - 2));
            if (dna_d[dna_index] == 'S')
                isEnd = 1;
            if (base2int(dna_d[dna_index]) != -1)
                cur_seq[k / 32] = (cur_seq[k / 32] << 2) % ((uint64) 1 << (2 * (k % 32))) +
                                  base2int(dna_d[dna_index]);
            if (isEnd)
                thread_dataR[list_index++] = UINT64_MAX;
            else
                thread_dataR[list_index++] = getHashValue(cur_seq, hash_b, k);
        }
    }
    delete [] cur_seq;
}

/*
 * Remove duplicate values from sorted hash value list.
*/
template<int BLOCK_THREADS, int ITEMS_PER_THREAD>
__device__ void BlockRemoveDupe (const int m, int *thread_dataS, int *thread_dataD, uint64 *thread_dataR,
                                int *thread_dataI, uint64 *sketch) {
    int offset = thread_dataS[0];
    for (int i = 0; i < ITEMS_PER_THREAD; ++i)
        if (thread_dataD[i] == 1 && offset < m)
            sketch[offset++] = thread_dataR[i];
}

/*
 * Get sketch back to input_d
*/
template<int BLOCK_THREADS, int ITEMS_PER_THREAD>
__device__ void BlockGetBack (const int m, uint64 *thread_dataR, int *thread_dataI, uint64 *sketch) {
    for (int i = 0; i < ITEMS_PER_THREAD; ++i)
        if (thread_dataI[i] < m)
            thread_dataR[i] = sketch[thread_dataI[i]];
}

/*
 * Core Function
 * Firstly, get hash value list from DNA sequence.
 *   For example, if dna_d is       [A,C,C,G,T,A,T,G,C,T,G,A,...].
 *   Then input_d should be         [6,5,6,1,0,2,4,2,1,2,6,3,...].
 * Secondly, block-sort the hash value list.
 *   For example, if input_d is     [6,5,6,1,0,2,4,2,1,2,6,3,...].
 *   Then output is still input_d:  [0,1,1,2,2,2,3,4,5,6,6,6,...].
 * Thirdly, mark the non-duplicate values.
 *   For example, if input_d is     [0,1,1,2,2,2,3,4,5,6,6,6,...].
 *   Then dupe_d should be          [1,1,0,1,0,0,1,1,1,1,0,0,...].
 * Then, scan dupe_d to find the offset of values in output.
 *   For example, if dupe_d is      [1,1,0,1,0,0,1,1,1,1,0,0,...].
 *   Then scan_d should be          [0,1,2,2,3,3,3,4,5,6,7,7,...].
 * Finally, get sketches of each block.
 */
template<int BLOCK_THREADS, int ITEMS_PER_THREAD>
__global__ void getBlockSketch (const int k, const int m, char *dna_d, uint64 *input_d,
                                int numElem_dna, int numElem_list, uint64 hash_b) {
    typedef cub::BlockStore<uint64, BLOCK_THREADS, ITEMS_PER_THREAD, cub::BLOCK_STORE_TRANSPOSE> BlockStoreR;
    typedef cub::BlockRadixSort <uint64, BLOCK_THREADS, ITEMS_PER_THREAD> BlockRadixSort;
    typedef cub::BlockDiscontinuity<uint64, BLOCK_THREADS> BlockDiscontinuity;
    typedef cub::BlockScan<int, BLOCK_THREADS> BlockScan;
    __shared__ union {
        typename BlockStoreR::TempStorage store;
        typename BlockRadixSort::TempStorage sort;
        typename BlockDiscontinuity::TempStorage dupe;
        typename BlockScan::TempStorage scan;
    } temp_storage;
    __shared__ uint64 sketch[BLOCK_THREADS * ITEMS_PER_THREAD];
    uint64 thread_dataR[ITEMS_PER_THREAD];
    int thread_dataD[ITEMS_PER_THREAD];
    int thread_dataS[ITEMS_PER_THREAD];
    int thread_dataI[ITEMS_PER_THREAD];

    int block_offset = blockIdx.x * (BLOCK_THREADS * ITEMS_PER_THREAD);
    int thread_offset = (BLOCK_THREADS * blockIdx.x + threadIdx.x) * ITEMS_PER_THREAD;

    // Get index for input list.
    for (int i = 0; i < ITEMS_PER_THREAD; ++i)
        thread_dataI[i] = 1;
    BlockScan(temp_storage.scan).ExclusiveSum(thread_dataI, thread_dataI);
    __syncthreads();

    BlockGetHashValues<BLOCK_THREADS, ITEMS_PER_THREAD>(k, dna_d, thread_offset, thread_dataR, hash_b);
    __syncthreads();
    BlockRadixSort(temp_storage.sort).Sort(thread_dataR);
    __syncthreads();
    BlockDiscontinuity(temp_storage.dupe).FlagHeads(thread_dataD, thread_dataR, cub::Inequality());
    __syncthreads();
    BlockScan(temp_storage.scan).ExclusiveSum(thread_dataD, thread_dataS);
    __syncthreads();
    BlockRemoveDupe<BLOCK_THREADS, ITEMS_PER_THREAD>(m, thread_dataS, thread_dataD, thread_dataR, thread_dataI, sketch);
    __syncthreads();
    BlockGetBack<BLOCK_THREADS, ITEMS_PER_THREAD>(m, thread_dataR, thread_dataI, sketch);
    __syncthreads();
    BlockStoreR(temp_storage.store).Store(input_d + block_offset, thread_dataR);
}

/* Merge.
 * Not Parallel Version Currently.
 */
void rMerge(const int m, uint64 *sketch_h, uint64 *output_h) {
    int pointer1 = 0, pointer2 = 0, count = 0;
    uint64 bucket[m];
    while (count < m) {
        if (sketch_h[pointer1] < output_h[pointer2]) {
            bucket[count++] = sketch_h[pointer1++];
        } else if (sketch_h[pointer1] > output_h[pointer2]) {
            bucket[count++] = output_h[pointer2++];
        } else if (sketch_h[pointer1] == output_h[pointer2]) {
            bucket[count++] = sketch_h[pointer1++];
            pointer2 ++;
        }
    }
    for (uint64 i = 0; i < m; i++)
        output_h[i] = bucket[i];
}

signature genSig(const int k, const int m, const int t, char *dnaList, int length, uint64 *hashes_b) {

    const int BLOCKS_NUM = 30;
    const int BLOCK_THREADS = 32 * 16;
    const int ITEMS_PER_THREAD = 4;

    // Compute CHUNKS_NUM and the start, end and record index.
    signature sig(t, vector<uint64>(m, UINT64_MAX));
    int CHUNKS_NUM;
    if (length % (BLOCKS_NUM * BLOCK_THREADS * ITEMS_PER_THREAD) == 0)
        CHUNKS_NUM = (length - k + 1) / (BLOCKS_NUM * BLOCK_THREADS * ITEMS_PER_THREAD);
    else
        CHUNKS_NUM = (length - k + 1) / (BLOCKS_NUM * BLOCK_THREADS * ITEMS_PER_THREAD) + 1;
    int *record = (int *) malloc(sizeof(int) * CHUNKS_NUM);
    int *start = (int *) malloc(sizeof(int) * CHUNKS_NUM);
    int *end = (int *) malloc(sizeof(int) * CHUNKS_NUM);
    for (int i = 0; i < (CHUNKS_NUM - 1); ++i) {
        record[i] = BLOCKS_NUM * BLOCK_THREADS * ITEMS_PER_THREAD + k - 1;
    }
    record[CHUNKS_NUM - 1] = length - (BLOCKS_NUM * BLOCK_THREADS * ITEMS_PER_THREAD) * (CHUNKS_NUM - 1);
    start[0] = 0;
    end[0] = record[0] - 1;
    if (CHUNKS_NUM >= 1)
        start[1] = record[0] - k + 1;
    for (int i = 1; i < CHUNKS_NUM - 1; i++) {
        end[i] = start[i] + record[i] - 1;
        start[i + 1] = end[i] + 1 - k + 1;
    }
    end[CHUNKS_NUM - 1] = length - 1;

    // Initialization.
    cudaError_t res;
    int numElem_dna = BLOCKS_NUM * BLOCK_THREADS * ITEMS_PER_THREAD + k - 1;
    int numElem_list = BLOCKS_NUM * BLOCK_THREADS * ITEMS_PER_THREAD;
    char *dna_h = (char *) malloc(sizeof(char) * numElem_dna);
    uint64 * output_h = (uint64 *) malloc(sizeof(uint64) * m);
    uint64 * sketch_h = (uint64 *) malloc(sizeof(uint64) * m);
    char *dna_d;
    uint64 * input_d;
    res = cudaMalloc(&dna_d, sizeof(char) * numElem_dna);
    CHECK(res);
    res = cudaMalloc(&input_d, sizeof(uint64) * numElem_list);
    CHECK(res);
    for (int i = 0; i < m; ++i)
        output_h[i] = UINT64_MAX;
    cout << "CHUNKS_NUM: " << CHUNKS_NUM << "  numElem_dna: " << numElem_dna << "  numElem_list: " << numElem_list << endl;

    for (int j = 0; j < t; j++) {
//        cout << "hash_index: " << j << "  hashes_b: " << hashes_b[j] << endl;
        for (int p = 0; p < CHUNKS_NUM; p++) {
//            cout << "\tchunk_index: " << p << endl;
            for (int i = 0; i < numElem_dna; i++) {
                if (i < record[p])
                    dna_h[i] = dnaList[i + start[p]];
                else
                    dna_h[i] = 'S';
            }
            res = cudaMemcpy((void *) (dna_d), (void *) (dna_h), numElem_dna * sizeof(char), cudaMemcpyHostToDevice);
            CHECK(res);

            getBlockSketch <BLOCK_THREADS, ITEMS_PER_THREAD> << < BLOCKS_NUM, BLOCK_THREADS >> >
                    (k, m, dna_d, input_d, numElem_dna, numElem_list, hashes_b[j]);

            res = cudaMemcpy((void *) (sketch_h), (void *) (input_d), m * sizeof(uint64), cudaMemcpyDeviceToHost);
            CHECK(res);
            rMerge(m, sketch_h, output_h);
        }
        for (int i = 0; i < m; i++) {
            sig[j][i] = output_h[i];
        }
    }
    cudaFree((void *) dna_d);
    cudaFree((void *) input_d);

    return sig;
}

