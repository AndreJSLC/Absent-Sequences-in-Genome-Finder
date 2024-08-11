#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <chrono>
#include <random>

using namespace std;
using namespace std::chrono;

// Function to read FASTA file
vector<pair<string, string>> read_fasta(const string& filepath) {
    vector<pair<string, string>> sequences;
    ifstream infile(filepath);
    string line, sequence, description;

    while (getline(infile, line)) {
        if (line[0] == '>') {
            if (!sequence.empty()) {
                transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
                sequences.push_back({ description, sequence });
                sequence.clear();
            }
            description = line.substr(1);
        }
        else {
            sequence += line;
        }
    }
    if (!sequence.empty()) {
        transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
        sequences.push_back({ description, sequence });
    }

    return sequences;
}
// Function to generate a random DNA sequence
string generate_random_sequence(int length, mt19937& rng) {
    string bases = "ACGT";
    string sequence;
    for (int i = 0; i < length; ++i) {
        sequence += bases[rng() % 4];
    }
    return sequence;
}

//// Function to calculate positional similarity using a sliding window
//bool calculate_similarity(const string& query, const string& target, int& max_score, double threshold) {
//    max_score = 0;
//    int query_len = query.size();
//    int target_len = target.size();
//
//    for (int i = 0; i <= target_len - query_len; ++i) {
//        int score = 0;
//        for (int j = 0; j < query_len; ++j) {
//            char query_base = query[j];
//            char target_base = target[i + j];
//            if (query_base == target_base ||
//                (target_base == 'R' && (query_base == 'A' || query_base == 'G')) ||
//                (target_base == 'Y' && (query_base == 'C' || query_base == 'T')) ||
//                (target_base == 'S' && (query_base == 'G' || query_base == 'C')) ||
//                (target_base == 'W' && (query_base == 'A' || query_base == 'T')) ||
//                (target_base == 'K' && (query_base == 'G' || query_base == 'T')) ||
//                (target_base == 'M' && (query_base == 'A' || query_base == 'C')) ||
//                (target_base == 'B' && (query_base == 'C' || query_base == 'G' || query_base == 'T')) ||
//                (target_base == 'D' && (query_base == 'A' || query_base == 'G' || query_base == 'T')) ||
//                (target_base == 'H' && (query_base == 'A' || query_base == 'C' || query_base == 'T')) ||
//                (target_base == 'V' && (query_base == 'A' || query_base == 'C' || query_base == 'G'))) {
//                score++;
//            }
//        }
//        if (static_cast<double>(score) / query_len * 100 > threshold) {
//            return true;  // Exit early if threshold is exceeded
//        }
//        max_score = max(max_score, score);
//    }
//
//    return false;  // Threshold not exceeded
//}
// Function to calculate positional similarity using a sliding window
bool calculate_similarity(const string& query, const string& target, int& max_score, double threshold) {
    int query_len = query.size();
    int target_len = target.size();

    // Calculate similarity for forward strand
    for (int i = 0; i <= target_len - query_len; ++i) {
        int score = 0;
        for (int j = 0; j < query_len; ++j) {
            if (query[j] == target[i + j]) {
                score++;
            }
        }
        if (static_cast<double>(score) / query_len * 100 > threshold) {
            return true;  // Exit early if threshold is exceeded
        }
        max_score = max(max_score, score);

    }

    // Calculate similarity for reverse complement strand
    string rev_comp_query = query;
    for (char& base : rev_comp_query) {
        switch (base) {
        case 'A': base = 'T'; break;
        case 'T': base = 'A'; break;
        case 'C': base = 'G'; break;
        case 'G': base = 'C'; break;
        }
    }
    reverse(rev_comp_query.begin(), rev_comp_query.end());

    for (int i = 0; i <= target_len - query_len; ++i) {
        int score = 0;
        for (int j = 0; j < query_len; ++j) {
            if (rev_comp_query[j] == target[i + j]) {
                score++;
            }
        }
        if (static_cast<double>(score) / query_len * 100 > threshold) {
            return true;  // Exit early if threshold is exceeded
        }
        max_score = max(max_score, score);
    }

    return false;  // Threshold not exceeded
}

int main() {
    string fastaFile;
    cout << "Enter the path to the genome FASTA file: ";
    getline(cin, fastaFile);



    // Start timing
    auto start = high_resolution_clock::now();
    vector<pair<string, string>> sequences = read_fasta(fastaFile);
    // Stop timing
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "Time taken to load genome into memory: " << duration.count() << " milliseconds" << endl;

    // Remove the key-value pairs where the key contains the substring "mitochondrion" because we don´t care about similarity in the mitochondrion
    sequences.erase(std::remove_if(sequences.begin(), sequences.end(),
        [](const pair<string, string>& p) {
            return p.first.find("mitochondrion") != string::npos;
        }),
        sequences.end());


    int numberOfSequencesToRetrieve;
    cout << "Enter the number of absent sequences you want the program to attempt to retrieve: ";
    cin >> numberOfSequencesToRetrieve;

    int sequenceLength;
    cout << "Enter the length of the sequences to generate(>12 nt): ";
    cin >> sequenceLength;

    int query_length = sequenceLength;

    unsigned int seed;
    cout << "Enter the seed for random sequence generation: ";
    cin >> seed;

    mt19937 rng(seed);

    // Get user-defined maximum homology threshold
    double max_homology_threshold;
    cout << "Enter the maximum homology percentage threshold: ";
    cin >> max_homology_threshold;

    // Ask the user if they want verbose output
    bool verbose;
    cout << "Do you want verbose output? (1 for yes, 0 for no): ";
    cin >> verbose;

    int absentSequencesFound = 0;
    ofstream outputFile("absent_sequences.txt");

    // Start timing
    start = high_resolution_clock::now();

    while (absentSequencesFound < numberOfSequencesToRetrieve) {
        string query_string = generate_random_sequence(sequenceLength, rng);
        bool homologyExceeded = false;
        cout << "Attempting sequence" + query_string << endl;

        for (const auto& seq : sequences) {
            int max_score = 0;
            if (calculate_similarity(query_string, seq.second, max_score, max_homology_threshold)) {
                homologyExceeded = true;
                if (verbose) {
                    cout << "Homology exceeded in " + seq.first << endl;
                }
                break;
            }
            if (verbose) {
                cout << "Max homology in " << seq.first << " = "
                    << (static_cast<double>(max_score) / query_length * 100)
                    << "%" << endl;
            }
        }

        if (!homologyExceeded) {
            absentSequencesFound++;
            outputFile << query_string << endl;
        }
    }

    outputFile.close();

    // Stop timing
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);
    cout << "Time taken to find absent sequences: " << duration.count() << " milliseconds" << endl;

    return 0;
}
