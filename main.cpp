#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <list>
#include <random>
#include <vector>
#include <string>

#include "Minhash.h"

using namespace std;

void output_signature(vector<vector<int>> sig1);

void help();

void usage();

int main(int argc, char *argv[]) {

    if (argc == 2 && string(argv[1]) == "help") help();
    if (argc < 4) usage();

    // Default values
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
    getline(file1, file_info1);
    getline(file2, file_info2);
    while (getline(file1, s1))
        if (!s1.empty())
            sequence1 += s1;
    while (getline(file2, s2))
        if (!s2.empty())
            sequence2 += s2;
    if (sequence1.size() < k || sequence2.size() < k) {
        cout << "k cannot be greater than the size of any document" << endl;
        exit(1);
    }
    file1.close();
    file2.close();

    clock_t ini_time;
    bool mode_found = false;
    double similarity, time;
    list<tuple<string, double, double>> results;
    if (cal_name == "all" || cal_name == "minhash") {
        if (t < 1) {
            cerr << endl;
            cerr << "You must provide a parameter --t=POSITIVE_INTEGER parameter for minhash modes!" << endl << endl;
            exit(1);
        }
        mode_found = true;
        ini_time = clock();
        vector<Hash> hashes = generateHashes(t, seed);
        vector<vector<int>> sig1 = generateSignature(k, m, sequence1, hashes);
        vector<vector<int>> sig2 = generateSignature(k, m, sequence2, hashes);
//        output_signature(sig1);
//        output_signature(sig2);
        similarity = computeSim(sig1, sig2);
        time = double(clock() - ini_time) / CLOCKS_PER_SEC;
        results.emplace_back("minhash", similarity, time);
    }
    if (!mode_found) usage();

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

void output_signature(vector<vector<int>> sig1) {
    for (int h = 0; h < sig1.size(); ++h) {
        cout << "sig1[" << h << "].size(): " << sig1[h].size() << "\t";
        for (int i = 0; i < sig1[h].size(); ++i) {
            cout << sig1[h][i] << " ";
        }
        cout << endl;
    }
}

void usage() {
    cout << "===========================" << endl;
    cerr << "Usage: " << endl << endl;
    cerr << "    ./SimCalculator FILE_ONE FILE_TWO MODE" << endl;
    cerr << endl;
    cerr << "    Possible MODEs are:" << endl;
    cerr << endl << bold_on;
    cerr << "        all" << endl;
    cerr << endl;
    cerr << "        jaccard" << endl;
    cerr << endl;
    cerr << "        jaccard_hash" << endl;
    cerr << endl;
    cerr << "        jaccard_roll" << endl;
    cerr << endl;
    cerr << "        minhash" << endl;
    cerr << endl;
    cerr << "        minhash_roll" << endl << bold_off;
    cout << endl;
    cerr << "Execute \"./SimCalculator help\" for an extended help section." << endl;
    cout << "===========================" << endl;
    exit(1);
}

void help() {
    cout << endl;
    cout << bold_on << "NAME" << bold_off << endl;
    cout << "    " << "SimCalculator" << endl;
    cout << endl;
    cout << bold_on << "USAGE" << bold_off << endl;
    cout << "    " << "SimCalculator FILE_ONE FILE_TWO " << bold_on << "MODE [PARAMETERS...]" << bold_off << endl;
    cout << endl;
    cout << "    " << "SimCalculator calculates the similarity between two text files FILE_ONE and FILE_TWO" << endl;
    cout << "    " << "and outputs it as a number between 0 and 1, where 1 means the two files are exactly" << endl;
    cout << "    " << "the same." << endl;
    cout << endl;
    cout << bold_on << "MODE" << bold_off << endl;
    cout << "    " << "There are five different modes which change the way SimCalculator computes the simi-" << endl;
    cout << "    " << "larity. Each may make use of different parameters, indicated as follows:" << endl;;
    cout << endl;
    cout << "    " << bold_on << "all" << bold_off << endl;
    cout << "        " << "This option executes all modes." << endl;
    cout << endl;
    cout << "    " << bold_on << "jaccard" << bold_off << endl;
    cout << "        " << "Performs the calculation by storing k-shingles as strings in a set for each file." << endl;
    cout << "        " << "Used parameters are: " << endl;
    cout << endl;
    cout << "            " << "--k=POSITIVE_INTEGER as shingle size" << endl;
    cout << endl;
    cout << "    " << bold_on << "jaccard_hash" << bold_off << endl;
    cout << "        " << "Performs the calculation by storing k-shingles as 4-byte hashes in a set for each" << endl;
    cout << "        " << "file. Uses less memory than jaccard, at the expense of a possible loss in" << endl;
    cout << "        " << "precision for large values of k. Used parameters are:" << endl;
    cout << endl;
    cout << "            " << "--k=POSITIVE_INTEGER as shingle size" << endl;
    cout << endl;
    cout << "    " << bold_on << "jaccard_roll" << bold_off << endl;
    cout << "        " << "Same as jaccard_hash, but hashes k-shingles by using a rolling hash. Better time" << endl;
    cout << "        " << "performance than jaccard_hash for greater values of k. Used parameters are:" << endl;
    cout << endl;
    cout << "            " << "--k=POSITIVE_INTEGER as shingle size" << endl;
    cout << endl;
    cout << "    " << bold_on << "minhash" << bold_off << endl;
    cout << "        " << "Calculates the similarity by computing minhash signatures for each document. Used" << endl;
    cout << "        " << "parameters are." << endl;
    cout << endl;
    cout << "            " << "--k=POSITIVE_INTEGER as shingle size" << endl;
    cout << endl;
    cout << "            " << "--t=POSITIVE_INTEGER" << bold_on << " (obligatory) " << bold_off
         << "as number of hash functions used" << endl;
    cout << endl;
    cout << "            " << "--seed=INTEGER as random generator seed" << endl;
    cout << endl;
    cout << "    " << bold_on << "minhash_roll" << bold_off << endl;
    cout << "        " << "Same as minhash mode, but uses a rolling hash. Better time performace tha minhash" << endl;
    cout << "        " << "for greater values of k. Used parameters are:" << endl;
    cout << endl;
    cout << "            " << "--k=POSITIVE_INTEGER as shingle size" << endl;
    cout << endl;
    cout << "            " << "--t=POSITIVE_INTEGER" << bold_on << " (obligatory) " << bold_off
         << "as number of hash functions used" << endl;
    cout << endl;
    cout << "            " << "--seed=INTEGER as random generator seed" << endl;
    cout << endl;
    cout << bold_on << "PARAMETERS" << bold_off << endl;
    cout << endl;
    cout << "    " << bold_on << "--k=POSITIVE_INTEGER" << bold_off << endl;
    cout << "        " << "Defaults to k=9. Indicates the size of the shingles used to calculate the simi-" << endl;
    cout << "        " << "larity between the documents." << endl;
    cout << endl;
    cout << "    " << bold_on << "--t=POSITIVE_INTEGER" << bold_off << endl;
    cout << "        " << "Indicates the number of hash functions used in minhash modes." << endl;
    cout << endl;
    cout << "    " << bold_on << "--seed=INTEGER" << bold_off << endl;
    cout << "        " << "Defaults to a random value. Used by minhash modes in their random generator number." << endl;
    cout << endl;
    cout << "    " << bold_on << "-e" << bold_off << endl;
    cout << "        " << "Output in experimentation format." << endl;
    cout << endl;
    exit(0);
}