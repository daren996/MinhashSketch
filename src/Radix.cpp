#include "Radix.h"
// 只用array
// merge的时候相邻merge
// √ 可以用8个桶，桶的大小也不必是uint64
// √ 提前检查重复的情况，这样就不必在merge的时候麻烦了
// √ 不要用pow，用<<更快

namespace rs {

    signature genSig_single(int k, int m, string &sequence, vector<Hash> &hashes) {
        int p = 1;
        signature sig(hashes.size(), vector<uint64>(m, UINT64_MAX));
        uint64 length = sequence.size();
        std::vector<uint64> record, start, end;
        for (int i = 0; i < (p - 1); ++i) {
            record.push_back(length / p);
        }
        record.push_back(length - length / p * (p - 1));
        start.push_back(0);
        end.push_back(record[0] - k);
        start.push_back(record[0] - k + 1);
        for (uint64 i = 1; i <= (p - 1); i++) {
            end.push_back(start[i] + record[i] - 1);
            if (i < (p - 1))
                start.push_back(start[i] + record[i + 1]);
        }
        for (int j = 0; j < hashes.size(); j++) {
            List list;
            for (uint64 i = 0; i < p; i++) {
                rSort(k, m, sequence, start[i], (end[i] + k - 1), hashes[j], list);
            }
            for (uint64 i = 0; i < m; i++) {
                sig[j][i] = list[i];
            }
        }
        return sig;
    }

    signature genSig_multi(int k, int m, string &sequence, vector<Hash> &hashes) {
        int p = 5;
        signature sig(hashes.size(), vector<uint64>(m, UINT64_MAX));

        // Compute the start and end index
        uint64 length = sequence.size();
        std::vector<uint64> record, start, end;
        for (int i = 0; i < (p - 1); ++i) {
            record.push_back(length / p);
        }
        record.push_back(length - length / p * (p - 1));
        start.push_back(0);
        end.push_back(record[0] - k);
        start.push_back(record[0] - k + 1);
        for (uint64 i = 1; i <= (p - 1); i++) {
            end.push_back(start[i] + record[i] - 1);
            if (i < (p - 1))
                start.push_back(start[i] + record[i + 1]);
        }

        // Start dispatching threads.
        for (int j = 0; j < hashes.size(); j++) {
            List list;
            std::vector<std::thread> tHolder;
            for (uint64 i = 0; i < p; i++) {
                rSort(k, m, sequence, start[i], (end[i] + k - 1), hashes[j], list);
//                std::thread t(rSort, k, m, std::ref(sequence), start[i], (end[i] + k - 1), std::ref(hashes[j]), std::ref(list));
//                tHolder.push_back(std::move(t));
            }
//            for (uint64 i = 0; i < tHolder.size(); ++i) {
//                tHolder[i].join();
//            }

            int blocks = p;
            while (blocks > 1) {
                std::vector<std::thread> tHolders;
                for (uint64 i = 0; i < (blocks / 2); i++) {
                    if (blocks % 2 == 0) {
                        std::thread t(rMerge, std::ref(list), start[i], start[i] + m - 1, start[i + blocks / 2],
                                      start[i + blocks / 2] + m - 1, m);
                        tHolders.push_back(std::move(t));
                    } else {
                        std::thread t(rMerge, std::ref(list), start[i], start[i] + m - 1, start[i + blocks / 2 + 1],
                                      start[i + blocks / 2 + 1] + m - 1, m);
                        tHolders.push_back(std::move(t));
                    }
                }
                for (uint64 i = 0; i < tHolders.size(); ++i) {
                    tHolders[i].join();
                }
                blocks = (blocks % 2 == 0) ? (blocks / 2) : (blocks / 2 + 1);
            }

            for (uint64 i = 0; i < m; i++) {
                sig[j][i] = list[i];
            }
        }
        return sig;
    }

    List rSort(int k, int m, string &sequence, uint64 begin, uint64 end, Hash &hash, List &list) {

        uint64 start = list.size();
        uint64 length = end - begin + 1;
        uint64 s_index = 0; // pointer of current base
        uint64 cur_seq[k / 32 + 1]; // current sub-sequence
        for (int i = 0; i < k / 32 + 1; ++i) {
            cur_seq[i] = 0;
        }

        // Get original list
        if (k < 32) {
            for (; s_index < k; ++s_index) {
                if (utils::base2int(sequence[s_index + begin]) != -1)
                    cur_seq[0] = (cur_seq[0] << 2) % ((uint64)1 << (2*k)) + utils::base2int(sequence[s_index + begin]);
                else
                    cerr << "ERROR:" << endl << "\t index: " << s_index + begin << endl << "\t base: "
                         << sequence[s_index + begin] << endl;
            }
            list.push_back(hash(cur_seq, (k / 32 + 1) * 8));
            for (; s_index < length; ++s_index) {
                if (utils::base2int(sequence[s_index + begin]) != -1)
                    cur_seq[0] = (cur_seq[0] << 2) % ((uint64)1 << (2*k)) + utils::base2int(sequence[s_index + begin]);
                else
                    cerr << "ERROR:" << endl << "\t index: " << s_index + begin << endl << "\t base: "
                         << sequence[s_index + begin] << endl;
                list.push_back(hash(cur_seq, (k / 32 + 1) * 8));
            }
        } else {
            for (; s_index < k; ++s_index) {
                if (utils::base2int(sequence[s_index + begin]) != -1)
                    cur_seq[s_index / 32] =
                            (cur_seq[s_index / 32] << 2) % UINT64_MAX + utils::base2int(sequence[s_index + begin]);
                else
                    cerr << "ERROR:" << endl << "\t index: " << s_index + begin << endl << "\t base: "
                         << sequence[s_index + begin] << endl;
            }
            list.push_back(hash(cur_seq, (k / 32 + 1) * 8));
            for (; s_index < length; ++s_index) {
                for (int i = 0; i < k / 32 - 1; ++i) {
                    cur_seq[i] = (cur_seq[i] << 2) + (cur_seq[i + 1] >> 62);
                }
                cur_seq[k / 32 - 1] = (cur_seq[k / 32 - 1] << 2) + (cur_seq[k / 32] >> ((k % 32) * 2 - 2));
                if (utils::base2int(sequence[s_index + begin]) != -1)
                    cur_seq[k / 32] = (cur_seq[k / 32] << 2) % ((uint64)1 << (2 * (k % 32))) +
                                      utils::base2int(sequence[s_index + begin]);
                else
                    cerr << "ERROR:" << endl << "\t index: " << s_index + begin << endl << "\t base: "
                         << sequence[s_index + begin] << endl;
                list.push_back(hash(cur_seq, (k / 32 + 1) * 8));
            }
        }

        rSortLowestToHighest(list, start, list.size() - 1);
        List temp;
        std::unordered_set<uint64> filter;
        for (int i = 0; i < m; ++i) {
            while (filter.insert(list[i + start]).second) {
                temp.push_back(list[i + start]);
            }
        }
        for (uint64 i = 0; i < m; ++i) {
            list[i + start] = temp[i];
        }
        return temp;
    }

    uint64 rMerge(List &list, uint64 begin1, uint64 end1, uint64 begin2, uint64 end2, int m) {
        uint64 pointer1 = begin1, pointer2 = begin2, count = 0;
        std::queue<uint64> bucket;
        while (pointer1 <= end1 && pointer2 <= end2 && count < m) {
            if (list[pointer1] < list[pointer2]) {
                bucket.push(list[pointer1]);
                count += 1;
                pointer1 += 1;
            } else if (list[pointer1] > list[pointer2]) {
                bucket.push(list[pointer2]);
                count += 1;
                pointer2 += 1;
            } else if (list[pointer1] == list[pointer2]) {
                bucket.push(list[pointer1]);
                count += 1;
                pointer1 += 1;
                pointer2 += 1;
            }
        }
        for (uint64 i = 0; i < count; i++) {
            list[i + begin1] = bucket.front();
            bucket.pop();
        }
        return count;
    }

    void rSortLowestToHighest(List &list, uint64 begin, uint64 end) {
        if (begin >= end) { return; }
        uint64 length = end - begin + 1;
        uint64 digits_max = 64;
        for (uint64 i = 0; i < digits_max; i+=3) {
            uint64 bucket[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            for (uint64 j = begin; j <= end; ++j) {
                uint64 curDigit = list[j] >> i & 7;
                bucket[curDigit + 1]++;
            }
            for (int j = 1; j < 9; j++) {
                bucket[j] += bucket[j - 1];
            }
            List tempL(length);
            for (uint64 j = begin; j <= end; ++j) {
                uint64 curDigit = list[j] >> i & 7;
                tempL[bucket[curDigit]++] = list[j];
            }
            for (uint64 j = 0; j < length; ++j) {
                list[j + begin] = tempL[j];
            }
        }
    }

}
