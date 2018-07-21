#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <list>
#include <random>
#include <vector>
#include <string>

#include "Radix.h"
#include "Help.h"

using namespace std;

void output_signature(vector<vector<uint64>> sig1);

// MinhashSketch.exe ../testing_files/sequence_clip1.fasta ../testing_files/sequence_clip2.fasta all -e --k=5 --m=10 --t=10
int main(int argc, char *argv[]) {

    if (argc == 2 && string(argv[1]) == "help") help();
    if (argc < 4) usage();

    // DEFAULT VALUES
    int k, m, t, seed;
    bool e;
    k = 9;
    m = 1; // the number of sketches
    t = 1; // the number of hash functions
    seed = random_device()();
    e = false;

    // PARSE FILE_ONE FILE_TWO MODE
    string name_one = string(argv[1]);
    string name_two = string(argv[2]);
    string cal_name = string(argv[3]);
    ifstream file1(name_one);
    if (file1.fail()) {
        std::cerr << "Unable to open file " << name_one << std::endl;
        exit(1);
    }
    ifstream file2(name_two);
    if (file2.fail()) {
        std::cerr << "Unable to open file " << name_two << std::endl;
        exit(1);
    }

    // PARSE PARAMETERS
    for (int i = 4; i < argc; ++i) {
        string param(argv[i]);
        if (param == "-e") {
            e = true;
        } else {
            int param_size = (uint) param.size();
            if (param_size >= 5) {
                auto index_eq = (uint) param.find('=');
                if (index_eq + 2 <= param_size) {
                    string param_name = param.substr(0, index_eq);
                    if (param_name == "--k") {
                        k = std::stoi(param.substr(index_eq + 1, param_size - index_eq - 1));
                    } else if (param_name == "--m") {
                        m = std::stoi(param.substr(index_eq + 1, param_size - index_eq - 1));
                    } else if (param_name == "--t") {
                        t = std::stoi(param.substr(index_eq + 1, param_size - index_eq - 1));
                    } else if (param_name == "--seed") {
                        seed = std::stoi(param.substr(index_eq + 1, param_size - index_eq - 1));
                    }
                }
            }
        }
    }
    if (k < 1) {
        std::cerr << "K value too small! Minimum: 1" << std::endl;
        exit(1);
    }
    if (m < 1) {
        std::cerr << "M value too small! Minimum: 1" << std::endl;
        exit(1);
    }

    // GET TWO SEQUENCES
    string file_info1, file_info2, sequence1, sequence2, s1, s2;
    utils::file_to_string(file1, file_info1, sequence1); // The first line is file information
    utils::file_to_string(file2, file_info2, sequence2);
    if (sequence1.size() < k || sequence2.size() < k) {
        cout << "k cannot be greater than the size of any document" << endl;
        exit(1);
    }
    file1.close();
    file2.close();

    // MAIN PROGRESS
    clock_t ini_time;
    bool mode_found = false;
    double similarity, time;
    list<tuple<string, double, double>> results;
    vector<Hash> hashes = generateHashes(t, seed);
    if (cal_name == "all" || cal_name == "minhash_regular") {
        if (t < 1) {
            cerr << endl;
            cerr << "You must provide a parameter --t=POSITIVE_INTEGER parameter for minhash modes!" << endl << endl;
            exit(1);
        }
        mode_found = true;
        ini_time = clock();
        vector<vector<uint64>> sig1 = rs::genSig_single(k, m, sequence1, hashes);
        vector<vector<uint64>> sig2 = rs::genSig_single(k, m, sequence2, hashes);
        cout << "sig1" << endl;
        output_signature(sig1);
        cout << "\nsig2" << endl;
        output_signature(sig2);
        cout << endl;
        similarity = computeSim(sig1, sig2);
        time = double(clock() - ini_time) / CLOCKS_PER_SEC;
        results.emplace_back("minhash_regular", similarity, time);
    }
    if (cal_name == "all" || cal_name == "minhash_parallel") {
        if (t < 1) {
            cerr << endl;
            cerr << "You must provide a parameter --t=POSITIVE_INTEGER parameter for minhash modes!" << endl << endl;
            exit(1);
        }
        mode_found = true;
        ini_time = clock();
        // vector<Hash> hashes = generateHashes(t, seed);
        vector<vector<uint64>> sig1 = rs::genSig_multi(k, m, sequence1, hashes);
        vector<vector<uint64>> sig2 = rs::genSig_multi(k, m, sequence2, hashes);
        cout << "sig1:  size:" << sig1[0].size() << endl;
        output_signature(sig1);
        cout << "\nsig2:  size:" << sig2[0].size() << endl;
        output_signature(sig2);
        cout << endl;
        similarity = computeSim(sig1, sig2);
        time = double(clock() - ini_time) / CLOCKS_PER_SEC;
        results.emplace_back("minhash_parallel", similarity, time);
    }
    if (!mode_found) usage();

    // OUTPUT RESULTS
    if (e) {
        cout << setw(12) << "cal_name" << setw(14) << "seed" << setw(5) << "k" << setw(5) << "m" << setw(7) << "t";
        cout << setw(13) << fixed << "time" << setw(13) << fixed << "similarity" << endl;
    } else {
        cout << "===========================" << endl;
        cout << "k:" << k << setw(7) << fixed << "m:" << m << setw(7) << fixed << "t:" << t << endl;
        cout << "===========================" << endl;
    }
    cout.precision(8);
    for (auto &result : results) {
        if (e) {
            cout << setw(12) << get<0>(result) << setw(14) << seed << setw(5) << k << setw(5) << m << setw(7) << t;
            cout << setw(13) << fixed << get<2>(result) << setw(13) << fixed << get<1>(result) << endl;
        } else {
            cout << uline_on << get<0>(result) << uline_off << endl;
            cout << "time: " << setw(21) << fixed << get<2>(result) << endl;
            cout << "similarity: " << setw(15) << fixed << get<1>(result) << endl;
            cout << "===========================" << endl;
        }
    }

}

void output_signature(vector<vector<uint64>> sig1) {
    for (int h = 0; h < sig1.size(); ++h) {
        cout << "sig[" << h << "].size(): " << sig1[h].size() << "\t";
        for (int i = 0; i < sig1[h].size(); ++i) {
            cout << hex << sig1[h][i] << dec << " ";
        }
        cout << endl;
    }
}
