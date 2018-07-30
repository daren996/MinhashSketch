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

using namespace std;

#define BLOCK_WIDTH 1 // 4 // 32
#define CEILING_DIVIDE(X, Y) (1 + (((X) - 1) / (Y)))
#define BLOCKS_NUM 5
#define CHECK(res) if(res!=cudaSuccess){printf("CHECK ERROR!\n");exit(-1);}

__global__
void rSort(int m, uint64 *input_d, uint64 *output_d, uint64 *record_d, uint64 *start_d, uint64 length) {

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
void rMerge(int m, uint64 *output_d, uint64 *record_d, uint64 *start_d, uint64 offset) {

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


signature genSig(int m, uint64 *list[], uint64 length, vector <Hash> &hashes) {

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
        res = cudaMemcpy((void *) (input_d), (void *) (input_h), length * sizeof(uint64), cudaMemcpyHostToDevice);CHECK(res);
        res = cudaMemcpy((void *) (output_d), (void *) (output_h), length * sizeof(uint64), cudaMemcpyHostToDevice);CHECK(res);
        res = cudaMemcpy((void *) (record_d), (void *) (record), BLOCKS_NUM * sizeof(uint64), cudaMemcpyHostToDevice);CHECK(res);
        res = cudaMemcpy((void *) (start_d), (void *) (start), BLOCKS_NUM * sizeof(uint64), cudaMemcpyHostToDevice);CHECK(res);

        rSort << < BLOCKS_NUM, 1 >> > (m, input_d, output_d, record_d, start_d, numElems);

        int offset = 2;
        while ((offset / 2) < BLOCKS_NUM) {
            rMerge << < BLOCKS_NUM, 1 >> > (m, output_d, record_d, start_d, offset);
            offset *= 2;
        }

        res = cudaMemcpy((void *) (output_h), (void *) (output_d), length * sizeof(uint64), cudaMemcpyDeviceToHost);CHECK(res);

        cudaFree((void *) input_d);
        cudaFree((void *) output_d);
        for (uint64 i = 0; i < m; i++) {
            sig[j][i] = output_h[i];
        }
    }
    return sig;

}


