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
//#define BLOCK_THREADS 16// 128
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

// Transfer DNA Sequence into Hash Value List.
template<int BLOCK_THREADS, int ITEMS_PER_THREAD>
__global__ void getList(const int k, char *dna_d, uint64 *input_d,
                        int numElem_dna, int numElem_list, uint64 hash_b) {
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int index = (BLOCK_THREADS * bx + tx) * ITEMS_PER_THREAD;
    if (index + k < numElem_dna) {
        int list_index = index;
        int dna_index = index;
        bool isEnd = 0;
        uint64 *cur_seq = new uint64[k / 32 + 1]; // current sub-sequence
        for (int i = 0; i < k / 32 + 1; ++i)
            cur_seq[i] = 0;
        if (k < 32) {
            for (; dna_index < index + k - 1; ++dna_index) {
                if (dna_d[dna_index] == 'S')
                    isEnd = 1;
                if (base2int(dna_d[dna_index]) != -1)
                    cur_seq[0] = (cur_seq[0] << 2) % ((uint64) 1 << (2 * k)) + base2int(dna_d[dna_index]);
            }
            for (; dna_index < index + ITEMS_PER_THREAD + k - 1; ++dna_index) {
                if (dna_d[dna_index] == 'S')
                    isEnd = 1;
                if (base2int(dna_d[dna_index]) != -1)
                    cur_seq[0] = (cur_seq[0] << 2) % ((uint64) 1 << (2 * k)) + base2int(dna_d[dna_index]);
                if (isEnd)
                    input_d[list_index++] = UINT64_MAX;
                else
                    input_d[list_index++] = getHashValue(cur_seq, hash_b, k);
            }
        } else {
            for (; dna_index < index + k; ++dna_index) {
                if (dna_d[dna_index] == 'S')
                    isEnd = 1;
                if (base2int(dna_d[dna_index]) != -1)
                    cur_seq[dna_index / 32] =
                            (cur_seq[dna_index / 32] << 2) % UINT64_MAX + base2int(dna_d[dna_index]);
            }
            if (isEnd)
                input_d[list_index++] = UINT64_MAX;
            else
                input_d[list_index++] = getHashValue(cur_seq, hash_b, k);
            for (; dna_index < index + ITEMS_PER_THREAD + k - 1; ++dna_index) {
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
                    input_d[list_index++] = UINT64_MAX;
                else
                    input_d[list_index++] = getHashValue(cur_seq, hash_b, k);
            }
        }
    }
}

// Radix Sort.
template<int BLOCK_THREADS, int ITEMS_PER_THREAD>
__global__ void rSort(uint64 *d_in) {
    // Specialize BlockLoad, BlockStore, and BlockRadixSort collective types
    typedef cub::BlockLoad<
            uint64, BLOCK_THREADS, ITEMS_PER_THREAD, cub::BLOCK_LOAD_TRANSPOSE> BlockLoadT;
    typedef cub::BlockStore<
            uint64, BLOCK_THREADS, ITEMS_PER_THREAD, cub::BLOCK_STORE_TRANSPOSE> BlockStoreT;
    typedef cub::BlockRadixSort <
            uint64, BLOCK_THREADS, ITEMS_PER_THREAD> BlockRadixSortT;
    // Allocate type-safe, repurposable shared memory for collectives
    __shared__ union {
        typename BlockLoadT::TempStorage load;
        typename BlockStoreT::TempStorage store;
        typename BlockRadixSortT::TempStorage sort;
    } temp_storage;
    // Obtain this block's segment of consecutive keys (blocked across threads)
    uint64 thread_keys[ITEMS_PER_THREAD];
    int block_offset = blockIdx.x * (BLOCK_THREADS * ITEMS_PER_THREAD);
    BlockLoadT(temp_storage.load).Load(d_in + block_offset, thread_keys);
    __syncthreads();    // Barrier for smem reuse
    // Collectively sort the keys
    BlockRadixSortT(temp_storage.sort).Sort(thread_keys);
    __syncthreads();    // Barrier for smem reuse
    // Store the sorted segment
    BlockStoreT(temp_storage.store).Store(d_in + block_offset, thread_keys);
}

/* Mark the Values Not Duplicate.
 * For example, if input_d is   [0,1,1,2,2,2,3,4,5,6,6,6,...].
 * Then dupe_d should be        [1,1,0,1,0,0,1,1,1,1,0,0,...].
 */
template<int BLOCK_THREADS, int ITEMS_PER_THREAD>
__global__ void getDupe (uint64 *input_d, int *dupe_d, int numElem_list) {
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int index = (BLOCK_THREADS * bx + tx) * ITEMS_PER_THREAD;
    if (index < numElem_list) {
        for (int i = 0; i < ITEMS_PER_THREAD; i++) {
            dupe_d[index + i] = 1;
            if (index != 0 && input_d[index + i] == input_d[index + i - 1])
                dupe_d[index + i] = 0;
        }
    }
}

/* Scan
 * For example, if dupe_d is    [1,1,0,1,0,0,1,1,1,1,0,0,...].
 * Then scan_d should be        [0,1,2,2,3,3,3,4,5,6,7,7,...].
 * Its duty is to find the offset of values in output.
 */
template<int BLOCK_THREADS, int ITEMS_PER_THREAD>
__global__ void getScan (int *scan_d, int *dupe_d, int numElem_list) {
    typedef cub::BlockLoad<
            int, BLOCK_THREADS, ITEMS_PER_THREAD, cub::BLOCK_LOAD_TRANSPOSE> BlockLoadT;
    typedef cub::BlockStore<
            int, BLOCK_THREADS, ITEMS_PER_THREAD, cub::BLOCK_STORE_TRANSPOSE> BlockStoreT;
    typedef cub::BlockScan<int, BLOCK_THREADS> BlockScanT;
    __shared__ union {
        typename BlockLoadT::TempStorage load;
        typename BlockStoreT::TempStorage store;
        typename BlockScanT::TempStorage scan;
    } temp_storage;
    int thread_data[4];
    int block_offset = blockIdx.x * (BLOCK_THREADS * ITEMS_PER_THREAD);
    BlockLoadT(temp_storage.load).Load(dupe_d + block_offset, thread_data);
    __syncthreads();
    BlockScanT(temp_storage.scan).ExclusiveSum(thread_data, thread_data);
    __syncthreads();
    BlockStoreT(temp_storage.store).Store(scan_d + block_offset, thread_data);
}


template<int BLOCK_THREADS, int ITEMS_PER_THREAD>
__global__ void getSketch (const int m, uint64 *input_d, int *scan_d, int *dupe_d,
                           uint64 *sketch_d, int numElem_list) {
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int index = (BLOCK_THREADS * bx + tx) * ITEMS_PER_THREAD;
    if (index < numElem_list) {
        int offset = scan_d[index];
        for (int i = 0; i < ITEMS_PER_THREAD; i++)
            if (dupe_d[index + i] == 1 && offset < m)
                sketch_d[offset++] = input_d[index + i];
    }
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

    const int BLOCKS_NUM = 5;
    const int BLOCK_THREADS = 16;// 128;
    const int ITEMS_PER_THREAD = 4;

    // Compute CHUNKS_NUM and the start, end and record index.
    signature sig(t, vector<uint64>(m, UINT64_MAX));
    int CHUNKS_NUM;
    if (length % (BLOCKS_NUM * BLOCK_THREADS * ITEMS_PER_THREAD) == 0)
        CHUNKS_NUM = (length - k + 1) / (BLOCKS_NUM * BLOCK_THREADS * ITEMS_PER_THREAD);
    else
        CHUNKS_NUM = (length - k + 1) / (BLOCKS_NUM * BLOCK_THREADS * ITEMS_PER_THREAD) + 1;
    //cout << "CHUNKS_NUM: " << CHUNKS_NUM << endl;
    //cout << "BLOCKS_NUM * BLOCK_THREADS * ITEMS_PER_THREAD: " << BLOCKS_NUM * BLOCK_THREADS * ITEMS_PER_THREAD << endl;
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
    //for (int i = 0; i < CHUNKS_NUM; i++) {
    //    cout << "start: " << start[i] << " | end: " << end[i] << " | record: " << record[i] << " " << endl;
    //}

    cudaError_t res;
    int numElem_dna = record[0];
    int numElem_list = BLOCKS_NUM * BLOCK_THREADS * ITEMS_PER_THREAD;
    char *dna_h = (char *) malloc(sizeof(char) * numElem_dna);
    uint64 * output_h = (uint64 *) malloc(sizeof(uint64) * m);
    uint64 * sketch_h = (uint64 *) malloc(sizeof(uint64) * m);
    char *dna_d;
    uint64 * input_d, *output_d, *sketch_d;
    int *dupe_d;
    int *scan_d;
    res = cudaMalloc(&dna_d, sizeof(char) * numElem_dna);
    CHECK(res);
    res = cudaMalloc(&input_d, sizeof(uint64) * numElem_list);
    CHECK(res);
    res = cudaMalloc(&dupe_d, sizeof(int) * numElem_list);
    CHECK(res);
    res = cudaMalloc(&scan_d, sizeof(int) * numElem_list);
    CHECK(res);
    res = cudaMalloc(&output_d, sizeof(uint64) * m);
    CHECK(res);
    res = cudaMalloc(&sketch_d, sizeof(uint64) * m);
    CHECK(res);
    for (int i = 0; i < m; i++)
        output_h[i] = UINT64_MAX;

    for (int j = 0; j < t; j++) {
        for (int p = 0; p < CHUNKS_NUM; p++) {
            for (int i = 0; i < numElem_dna; i++) {
                if (i < record[p])
                    dna_h[i] = dnaList[i + start[p]];
                else
                    dna_h[i] = 'S';
            }
            res = cudaMemcpy((void *) (dna_d), (void *) (dna_h), numElem_dna * sizeof(char), cudaMemcpyHostToDevice);
            CHECK(res);
//            res = cudaMemcpy((void *) (output_d), (void *) (output_h), m * sizeof(uint64), cudaMemcpyHostToDevice);
//            CHECK(res);

            getList<BLOCK_THREADS, ITEMS_PER_THREAD> << < BLOCKS_NUM, BLOCK_THREADS >> >
                    (k, dna_d, input_d, numElem_dna, numElem_list, hashes_b[j]);
            rSort<BLOCK_THREADS, ITEMS_PER_THREAD>  << < BLOCKS_NUM, BLOCK_THREADS >> >
                    (input_d);
            getDupe <BLOCK_THREADS, ITEMS_PER_THREAD> << < BLOCKS_NUM, BLOCK_THREADS >> >
                    (input_d, dupe_d, numElem_list);
            getScan <BLOCK_THREADS, ITEMS_PER_THREAD> << < BLOCKS_NUM, BLOCK_THREADS >> >
                    (scan_d, dupe_d, numElem_list);
            getSketch <BLOCK_THREADS, ITEMS_PER_THREAD> << < BLOCKS_NUM, BLOCK_THREADS >> >
                    (m, input_d, scan_d, dupe_d, sketch_d, numElem_list);

            res = cudaMemcpy((void *) (sketch_h), (void *) (sketch_d), m * sizeof(uint64), cudaMemcpyDeviceToHost);
            CHECK(res);
            rMerge(m, sketch_h, output_h);
//            res = cudaMemcpy((void *) (output_h), (void *) (sketch_d), m * sizeof(uint64), cudaMemcpyDeviceToHost);
//            CHECK(res);
        }
        for (int i = 0; i < m; i++) {
            sig[j][i] = output_h[i];
        }
    }
    cudaFree((void *) dna_d);
    cudaFree((void *) input_d);
    cudaFree((void *) output_d);
    cudaFree((void *) dupe_d);
    cudaFree((void *) scan_d);
    cudaFree((void *) sketch_d);

    return sig;
}




















/*
 * The following is old program accepting uint64* input.
 * */
#define BLOCKS_NUM 5
__global__
void rSort_list(int m, uint64 *input_d, uint64 *output_d, uint64 *record_d, uint64 *start_d, uint64 length) {

    int block_i = blockIdx.x;

    if (block_i < BLOCKS_NUM) {
        uint64 digits_max = 64;
        for (uint64 i = 0; i < digits_max; i += 3) {
            uint64 bucket[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            for (uint64 j = start_d[block_i]; j < start_d[block_i] + record_d[block_i]; ++j) {
                uint64 curDigit = input_d[j] >> i & 7;
                bucket[curDigit + 1]++;
            }
            for (int j = 1; j < 9; j++) {
                bucket[j] += bucket[j - 1];
            }
            for (uint64 j = start_d[block_i]; j < start_d[block_i] + record_d[block_i]; ++j) {
                uint64 curDigit = input_d[j] >> i & 7;
                output_d[start_d[block_i] + bucket[curDigit]++] = input_d[j];
            }
            for (uint64 j = start_d[block_i]; j < start_d[block_i] + record_d[block_i]; ++j) {
                input_d[j] = output_d[j];
            }
        }
        // Remove duplicate
        for (uint64 i = 0, j = start_d[block_i] + 1; j < start_d[block_i] + record_d[block_i] && i < m; ++j)
            if (output_d[start_d[block_i] + i] != input_d[j])
                output_d[++i + start_d[block_i]] = input_d[j];
    }

}

__global__
void rMerge_list(int m, uint64 *output_d, uint64 *record_d, uint64 *start_d, uint64 offset) {

    int block_i = blockIdx.x;

    if (block_i + (offset / 2) < BLOCKS_NUM && block_i % offset == 0) {
        uint64 begin1, begin2, end1, end2;
        begin1 = start_d[block_i];
        end1 = start_d[block_i] + record_d[block_i] - 1;
        begin2 = start_d[block_i + (offset / 2)];
        end2 = start_d[block_i + (offset / 2)] + record_d[block_i + (offset / 2)] - 1;
        uint64 pointer1 = begin1, pointer2 = begin2, count = 0;
        uint64 *bucket = (uint64 *) malloc(m * sizeof(uint64));
        while (pointer1 <= end1 && pointer2 <= end2 && count < m) {
            if (output_d[pointer1] < output_d[pointer2]) {
                bucket[count++] = output_d[pointer1];
                pointer1 += 1;
            } else if (output_d[pointer1] > output_d[pointer2]) {
                bucket[count++] = output_d[pointer2];
                pointer2 += 1;
            } else if (output_d[pointer1] == output_d[pointer2]) {
                bucket[count++] = output_d[pointer1];
                pointer1 += 1;
                pointer2 += 1;
            }
        }
        for (uint64 i = 0; i < count; i++)
            output_d[i + begin1] = bucket[i];
    }

}

signature genSig_list(int m, uint64 *list[], uint64 length, vector <Hash> &hashes) {

    signature sig(hashes.size(), vector<uint64>(m, UINT64_MAX));

    // Compute the start and end index
    uint64 *record = (uint64 *) malloc(sizeof(uint64) * BLOCKS_NUM);
    uint64 *start = (uint64 *) malloc(sizeof(uint64) * BLOCKS_NUM);
    for (int i = 0; i < (BLOCKS_NUM - 1); ++i) {
        record[i] = length / BLOCKS_NUM + 1;
    }
    record[BLOCKS_NUM - 1] = length - (length / BLOCKS_NUM + 1) * (BLOCKS_NUM - 1);
    start[0] = 0;
    for (uint64 i = 1; i < BLOCKS_NUM; i++) {
        start[i] = start[i - 1] + record[i - 1];
    }

    // Start dispatching blocks
    for (int j = 0; j < hashes.size(); j++) {
        cudaError_t res;
        uint64 numElems = length / BLOCKS_NUM + 1;
        uint64 *input_d, *output_d, *record_d, *start_d;
        uint64 *input_h = (uint64 *) malloc(sizeof(uint64) * length);
        uint64 *output_h = (uint64 *) malloc(sizeof(uint64) * length);
        res = cudaMalloc(&input_d, sizeof(uint64) * length);
        CHECK(res);
        res = cudaMalloc(&output_d, sizeof(uint64) * length);
        CHECK(res);
        res = cudaMalloc(&record_d, sizeof(uint64) * BLOCKS_NUM);
        CHECK(res);
        res = cudaMalloc(&start_d, sizeof(uint64) * BLOCKS_NUM);
        CHECK(res);
        for (uint64 i = 0; i < length; i++) {
            input_h[i] = list[j][i];
            output_h[i] = UINT64_MAX;
        }
        res = cudaMemcpy((void *) (input_d), (void *) (input_h), length * sizeof(uint64), cudaMemcpyHostToDevice);
        CHECK(res);
        res = cudaMemcpy((void *) (output_d), (void *) (output_h), length * sizeof(uint64), cudaMemcpyHostToDevice);
        CHECK(res);
        res = cudaMemcpy((void *) (record_d), (void *) (record), BLOCKS_NUM * sizeof(uint64), cudaMemcpyHostToDevice);
        CHECK(res);
        res = cudaMemcpy((void *) (start_d), (void *) (start), BLOCKS_NUM * sizeof(uint64), cudaMemcpyHostToDevice);
        CHECK(res);

        rSort_list << < BLOCKS_NUM, 1 >> > (m, input_d, output_d, record_d, start_d, numElems);

        int offset = 2;
        while ((offset / 2) < BLOCKS_NUM) {
            rMerge_list << < BLOCKS_NUM, 1 >> > (m, output_d, record_d, start_d, offset);
            offset *= 2;
        }

        res = cudaMemcpy((void *) (output_h), (void *) (output_d), length * sizeof(uint64), cudaMemcpyDeviceToHost);
        CHECK(res);

        cudaFree((void *) input_d);
        cudaFree((void *) output_d);
        for (uint64 i = 0; i < m; i++) {
            sig[j][i] = output_h[i];
        }
    }
    return sig;

}




