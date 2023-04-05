#include <seqan/index.h>
#include <seqan/seq_io.h>

//gm
#include <experimental/filesystem>
#include <vector>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <math.h>
#include <cstdlib>
#include <typeinfo>
//#include "omp.h"
//gm

using namespace std;
using namespace seqan;

typedef Iterator<StringSet<DnaString> >::Type TStringSetIterator;


// Initializing globals:
int exactMatches = 0, eta_seconds = 0, noOfThreads;
int done_datasets=0, done_trials=0, done_buckets=0;
int l, d, k, s, m, datasets, seqidx, tmp_round, tmp_trial, buckets_quantity;
double avrg_performance_coefficient = 0;

chrono::steady_clock::time_point program_start, datasets_start, step;

vector<int> time_left(3,0);
vector<double> times;
vector<thread> threads;
vector<vector<string>> parameters, pm;
mutex mtx_done_buckets, mtx_bucket_conseqs;

map<int, vector<pair<int,DnaString>>> foundMatches;
StringSet<DnaString> sequences;

bool printProgress = false;


//gm

/**
 * GENMAP FUNCTION
 *
 * GenMap function to load raw files
 * https://github.com/cpockrandt/genmap/wiki/#how-to-load-raw-files-map-freq8-freq16-in-c
 *
 * @tparam value_t
 * @param vec
 * @param path
 */
template <typename value_t>
void load(vector<value_t> & vec, string path)
{
    ifstream file(path, ios::binary);
    if (!file.eof() && !file.fail())
    {
        file.seekg(0, ios_base::end);
        streampos fileSize = file.tellg();
        vec.resize(fileSize / sizeof(value_t));
        file.seekg(0, ios_base::beg);
        file.read(reinterpret_cast<char*>(& vec[0]), fileSize);
        file.close();
        return;
    }
    // something went wrong ...
}
/**
 * GENMAP FUNCTION
 *
 * @param path_filename
 * @param filename
 * @param motif_length
 * @param mismatches
 * @return
 */
vector<uint8_t> getGenMapFrequencyVector(string path_filename, string filename, int motif_length, int mismatches) {
    //mismatches = 2;
    //$ ./genmap index -F /path/to/fasta.fasta -I /path/to/index/folder
    vector<uint8_t> frequency_vector;
    string output_folder_name = "_output_";
    string genmap_command =
            "genmap index -F " + path_filename+".fasta" + " -I ./_indeces_;"
            + " genmap map "
            + " -E " + to_string(mismatches)
            + " -K " + to_string(motif_length) + " -I ./_indeces_"
            + " -O ./_output_ -r -fs";
    //mkdir(output_folder_name.c_str(), 777);
    size_t ftpos = filename.find(".fasta");
    filename = filename.substr(0, ftpos);
    string mkdir = "mkdir _output_";
    system(mkdir.c_str());
    cout << "running " << genmap_command << endl;
    system(genmap_command.c_str());
    load(frequency_vector, "./_output_/" + filename + ".genmap.freq8");
    //experimental::filesystem::remove_all("./_output_/");
    //experimental::filesystem::remove_all("./_indeces_/");
    system("rm -rf ./_output_");
    system("rm -rf ./_indeces_");
    cout << "length of 1. fv: " << length(frequency_vector) << endl;
    cout << "./_output_/" << filename << ".freq8" << endl;
    return frequency_vector;
}
vector<uint8_t> getGenMapFrequencyVectorOPS(string path_to_directory, string filename, int motif_length, int mismatches) {
    //mismatches = 3;
    //"genmap index -FD /path/to/directory/with/fasta/files -I /multi/fasta/index/output"
    //"genmap map -I /multi/fasta/index/output"
    vector<uint8_t> frequency_vector;
    string output_folder_name = "_output_";
    string genmap_command =
            "genmap index -FD " + path_to_directory + " -I ./_indeces_;"
            + " genmap map "
            + " -E " + to_string(mismatches)
            + " -K " + to_string(motif_length) + " -I ./_indeces_"
            + " -O ./_output_ -r -fs -ep";
    //mkdir(output_folder_name.c_str(), 777);
    //size_t pos = filename.find("genmap_");
    size_t pos = length("genmap_");
    size_t len = length(filename);
    filename = filename.substr(pos, len);
    //cout << "then it is : " << filename << endl;
    string mkdir = "mkdir _output_";
    system(mkdir.c_str());
    //cout << "running " << genmap_command << endl;
    system(genmap_command.c_str());
    load(frequency_vector, "./_output_/" + filename + ".genmap.freq8");
    //experimental::filesystem::remove_all("./_output_/");
    //experimental::filesystem::remove_all("./_indeces_/");
    cout << "length of 1. fv: " << length(frequency_vector) << endl;
    cout << "./_output_/" << filename << ".freq8" << endl;
    system("rm -rf ./_output_");
    system("rm -rf ./_indeces_");
    return frequency_vector;
}
/**
 * HELPER FUNCTION
 *
 * Remove duplicates
 * https://www.techiedelight.com/remove-duplicates-vector-cpp/
 *
 * @param v
 * @return
 */
vector<int> removeDuplicates(vector<int> v) {
    sort(v.begin(), v.end());
    v.erase(unique(v.begin(), v.end()), v.end());
    return v;
}

//gm





/**
 * The maps in this program have a vector at each key. Adding a new element to key x means
 * "add element to the vector at key x" or, if key x has not been initialized yet, "add a
 * vector to key x, then add the element to that vector". This helper function allows us
 * to define the type of the element to be added.
 */
template<typename T> void addToMap(map<int, vector<T>>& map, int key, T item) {
    if(map.count(key) > 0) //if it already exist
        map[key].push_back(item);
    else
        map.insert(pair<int, vector<T>>(key,vector<T> (1,item)));
}


//gm

/**
 * HELPER FUNCTION
 *
 *
 * @param candidates
 */

vector<int> covertFreqVecToIntVec(vector<uint8_t> frequency_vector){
    vector<int> frequency_vector_int;
    for (int i = 0; i < length(frequency_vector); i++) {
        frequency_vector_int.push_back((int) frequency_vector[i]);
        //cout << (int) frequency_vector[i] << " ";
        if((int) frequency_vector[i] == 20) cout << "20!!!!!!!!!!!!!" << endl;
    }
    return frequency_vector_int;
}

//gm




/**
 * A time vector in this program is a vector with 3 ints:
 * Index 0 is the hours, 1 the minutes, 2 the seconds.
 * This helper function returns a string in the format of hh:mm:ss.
 */
string timeVectorToString(vector<int>t) {
    string result;
    for(int i=0; i<t.size(); i++)
        result+= ((t[i]<10) ? "0":"") + to_string(t[i]) + ((i<t.size()-1) ? ":":"");
    return result;
}



/**
 * This helper function parses a csv file into a 2 dimensional vector of strings.
 * The delimiter can be specified optionally, the default is the tab character.
 * You can also specify from which line on the csv should be parsed, as many csv
 * files have headers on line 0 (eg the parameters csv for this program).
 * Cell at row1 col5 would be accessed with vector[0][4]
 */
vector<vector<string>> csvParseIntoVector(istream& stream, char delim = '\t', int startRow = 0) {
    vector<vector<string>> result;
    string line;
    for(int i=0; i<startRow; i++) // Skip lines in csv
        getline(stream,line);
    while(getline(stream,line)) { // parse every other line
        vector<string> parsedRow;
        stringstream row(line);
        string cell;
        while(getline(row, cell, delim))
            parsedRow.push_back(cell);
        if (!row && cell.empty())
            parsedRow.push_back(""); // add an empty string if trailing delim
        result.push_back(parsedRow);
    }
    return result;
}



/**
 * This function iterates through a list of consensus sequences
 * and returns the consensus sequence with the lowest score.
 * If a score of 0 is found, the consensus sequence is returned right away,
 * if not, the first consensus sequence with the lowest score available score
 * is returned.
 */
pair<DnaString,int> best_consensus_of(vector<pair<DnaString,int>>& conseqs) {
    pair<DnaString,int> best_conseq;
    int cur_min_score = -1;
    for(auto conseq : conseqs) {
        if(conseq.second == 0) // lowest possible score, return right away
            return conseq;
        if(cur_min_score == -1 || conseq.second < cur_min_score) { 	// for the first run, set cur_min_score to first conseq score
            cur_min_score = conseq.second; // from the first run take the smallest
            best_conseq = conseq;
        }
    }
    return best_conseq;
}



/**
 * datasets a double x to n decimals.
 */
double rwp(double num, int precision) { //round with precision = rwp
    double factor = pow(10.0,(double) precision+1);
    return (int)(((num * factor)+5)/10)/(factor/10);
}






/**
 * Get the hamming distance of two DnaStrings.
 * from https://www.geeksforgeeks.org/hamming-distance-two-strings/
 */
int hammingDist(DnaString& str1, DnaString& str2) {
    int i = 0, count = 0;
    for(int x = 0; x < length(str1); x++) {
        if (str1[i] != str2[i])
            count++;
        i++;
    }
    return count;
}



/**
 * Return a substr of a DnaString.
 */
DnaString substr(DnaString& seq, int start, int length) {
    DnaString res;
    for(int i=start; i<start+length; i++)
        res += seq[i];
    return res;
}









/**
 * PROJECTION FUNCTION
 *
 *
 *
 */
void addBackgroundProbability(float* background_probability_distribution, vector<vector<float>>& Wh, int& motif_length, int number_of_sequences) {
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < motif_length; j++)
            Wh[i][j] = Wh[i][j] / number_of_sequences + background_probability_distribution[i]; //in percent + Laplace correction
}

/**
 * Initializing weight matrix
 */
void initWh(int& motif_length, pair<int,vector<pair<int,int>>>& bucket, vector<vector<float>>& Wh) {
    Wh = vector<vector<float>>(4,vector<float>(motif_length,0));
    float P[4] = {0.25, 0.25, 0.25, 0.25}; //background probability distribution

    for(auto p : bucket.second) { //for each element in the bucket
        int i = p.first, j = p.second; // i = seqNr , j = position of l-mer in sequence
        for (int u = 0; u < motif_length; u++) //calculate the frequencies of the bases at every position u in the l-mers
            Wh[ordValue(sequences[i][u+j])][u]++; //ordValue turns a nt into an int
    }

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < motif_length; j++)
            Wh[i][j] = Wh[i][j]/bucket.second.size() + P[i]; //in percent + Laplace correction
}




void refine(int& motif_length, vector<vector<float>>& Wh, vector<vector<float>>& posM) {
    vector<vector<float>> ref_Wh_tmp(Wh);
    vector<float> sums = vector<float>(motif_length,0);

    for (TStringSetIterator it = begin(sequences); it != end(sequences); ++it){
        int seq = distance(begin(sequences),it); //seqNr
        float denominator = 0; //denominator
        for (int i = 0; i < length(sequences[0]) - motif_length + 1; ++i){
            float numerator = 1;
            for (int u = 0; u < motif_length; u++)
                numerator *= Wh[ordValue((*it)[u+i])][u];
            posM[seq][i] = numerator; //numerator
            denominator += numerator;
        }
        for(int i = 0; i < length(sequences[0]) - motif_length + 1; ++i)
            posM[seq][i] /= denominator;
    }

    for(int start = 0; start < motif_length; start++)
        for(int i=0; i < length(sequences); i++)
            for (int j=start, pos=0; (pos < length(sequences[0]) - motif_length + 1) ; j++, pos++)
                ref_Wh_tmp[ordValue(sequences[i][j])][start] += posM[i][pos]; //numerator
    for(int i = 0; i<4; i++) //iterating through Wh
        for(int j = 0; j < motif_length; j++)
            sums[j] += ref_Wh_tmp[i][j]; //for the denominator
    for(int i = 0; i<4; i++)
        for(int j = 0; j < motif_length; j++)
            Wh[i][j] = ref_Wh_tmp[i][j] / sums[j];
}


/**
 * create the stringSet T of the l-mers
 * take an l-mer from each sequence using the posM
 * (fill T with this bucket's lmers)
 */
void get_consensus_seq(int& motif_length, int& d,vector<vector<float>>& posM, vector<pair<DnaString,int>>& bucket_conseqs) {

    StringSet<DnaString> T;
    DnaString conseq_of_t; //consensus sequnce
    int score_of_T = 0; // the number of l-mers whose hamming distance to the consensus sequence is larger than d

    for(int i = 0, row = 0, col = 0; i < length(sequences); i++) {
        DnaString current_lmer;
        float local_max = 0;
        for(int j = 0; j < length(sequences[0]) - motif_length + 1; j++) {
            if(posM[i][j] > local_max) {
                local_max = posM[i][j];
                row = i;
                col = j;
            }
        }
        current_lmer = substr(sequences[row],col,motif_length);
        appendValue(T, current_lmer);
    }

    for(int j=0; j<motif_length;j++) { ;
        int scores[4] = {0};
        int max_score = 0;
        char character;
        for(int i=0; i<length(T); i++) {
            scores[ordValue(T[i][j])]++;
            if(max_score < scores[ordValue(T[i][j])]) {
                max_score = scores[ordValue(T[i][j])];
                character = T[i][j];
            }
        }
        conseq_of_t += character;
    }

    for(auto lmer : T)
        if(hammingDist(conseq_of_t, lmer) > d)
            score_of_T++;

    mtx_bucket_conseqs.lock();
    bucket_conseqs.push_back(pair<DnaString,int>(conseq_of_t,score_of_T));
    mtx_bucket_conseqs.unlock();
}

//gm
map<int, vector<pair<int,int>>> convertToMap(vector<pair<int,int>>& vector_of_pairs) {
    map<int, vector<pair<int,int>>> vector_of_pairs_as_map;
    for(int i = 0; i < vector_of_pairs.size(); i++) {
        addToMap<pair<int,int>>(vector_of_pairs_as_map, vector_of_pairs.at(i).first, vector_of_pairs.at(i));
    }
    return vector_of_pairs_as_map;
}
map<int, vector<pair<int,int>>> convertPairToMap(pair<int,int> pair_to_convert) {
    map<int, vector<pair<int,int>>> pair_as_map;
    addToMap<pair<int,int>>(pair_as_map, pair_to_convert.first, pair_to_convert);
    return pair_as_map;
}
vector<pair<int,int>> convertFromMap(map<int, vector<pair<int, DnaString>>>& vector_of_pairs_as_map) {
    vector<pair<int,int>> vector_of_pairs;
    for(pair<int, vector<pair<int,DnaString>>> row : vector_of_pairs_as_map){
        int sequence_number = row.first;
        for(int i = 0; i < row.second.size(); i++) {
            vector_of_pairs.push_back(pair<int,int>((int) row.first, (int) row.second.at(i).first));
        }
    }
    return vector_of_pairs;
}
//gm



/**
 * Convert an int of seconds into a vector of three ints.
 * Index 0 for hours, 1 minutes and 2 seconds.
 */
vector<int> secondsToHours(int seconds) {
    int hours = seconds / 3600; seconds -= (3600*hours);
    int minutes = seconds / 60; seconds -= (60*minutes);
    return vector<int> {hours,minutes,seconds};
}



/**
 * Print the progress of the process.
 */
void print_progress() {
    while(printProgress) {
        this_thread::sleep_for(chrono::milliseconds(200)); // refresh rate: around 5 times per second
        double trial_progress = 100*(double) done_buckets / (double) buckets_quantity;
        double dataset_progress = (100*(double) (done_trials%m) / (double) m)+(trial_progress/m);
        double overall_progress = (100*(double) done_datasets/(double) datasets)+(dataset_progress/datasets);
        int curdec = (int) rwp(trial_progress,0);
        int diff = chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-datasets_start).count(); //now-program_start
        int total_seconds = (100.0 / overall_progress)*diff;
        time_left = secondsToHours(total_seconds - diff);

        cout << "\rOAP: " << rwp(overall_progress,0) << "%    "; // Overall Progress
        cout << "ETA(hh:mm:ss): " << timeVectorToString(time_left) << "    "; // Estimated time left
        cout << "Dataset [" << tmp_round <<"/"<<datasets << "]: " << rwp(dataset_progress,0) << "%    "; // Dataset progress
        cout << "Trial ["<<tmp_trial<<"/"<<m<<"]: "<< curdec <<"% ["; // Trial progress
        for(int i=0; i<curdec/5;i++) cout << "#"; // Trial progress bar, reached positions filled with #
        for(int i= curdec/5; i<20; i++) cout << " "; // Trial progress bar, left positions filled with spaces
        cout << "]       "<< flush; // add some spaces to full override previous line, flush the stream onto the console
    }
}



/**
 * For each bucket that has a size euqals or greater than the threshold,
 * compute and refine weight and position matrix. Then save the consesus sequence
 * of the bucket. This function is designed to be run by multiple threads. Each thread
 * gets a range of buckets to process. That range is bstart till bend.
 */
bool motif_refinement(int motif_length, int threshold, int trial, int bstart, int bend, map<int, vector<pair<int,int>>>& buckets, vector<pair<DnaString,int>>& bucket_conseqs) {
    int i = 0;
    for(pair<int,vector<pair<int,int>>> bucket : buckets){
        if(i>=bstart && i<=bend) {
            if(bucket.second.size() >= threshold){ //number of pairs = number of elements in the bucket
                vector<vector<float>> Wh;
                vector<vector<float>> posM = vector<vector<float>>(length(sequences),(vector<float>(length(sequences[0])-motif_length+1,0)));
                initWh(motif_length, bucket, Wh); //initialize weight matrix Wh for prob of a base in the motif
                for(int refine_iter = 0; refine_iter < 5; refine_iter++) //Refine weight matrix W and posM until convergence
                    refine(motif_length, Wh, posM);
                get_consensus_seq(motif_length, d, posM, bucket_conseqs); //CONSENSUS SEQUENCE
            }
            mtx_done_buckets.lock(); done_buckets++; mtx_done_buckets.unlock();
        }
        i++;
    }
    return true;
}




//gm
double calculatePerformanceCoefficient(DnaString pattern, double lmer_performance_coefficient, int exact_matches) {
    int z = 0;
    int n = 0;
    if(DnaString(pm[0][2]) == pattern) {
        z = 1; n = 1; exact_matches++;
    } else {
        for(int i = 0; i < length(sequences); i++) {
            int planted_position = stoi(pm[i+1][1]);
            if(foundMatches.count(i) > 0) {
                int lowerBound = planted_position - (l - 1);
                int upperBound = planted_position + (l - 1);
                int u = l;
                int current_fm;
                for(pair<int, DnaString> foundMatch: foundMatches[i]) {
                    current_fm = foundMatch.first;
                    if(current_fm >= lowerBound && current_fm <= upperBound)
                        break;
                }
                int tail = (current_fm >= lowerBound && current_fm <= upperBound) ? abs(planted_position - current_fm) : l;
                u = min(u, tail);
                z += (l - u);
                n += (l + u);
            } else {
                n += l;
            }
        }
    }
    lmer_performance_coefficient += (double)z / (double)n;
    //cout << "The performance coefficient is " << roundWithPrecision((double)z / (double)n, 2) << "." << endl;
    return lmer_performance_coefficient;
}

DnaString runGenMap2(int motif_length, vector<pair<int,int>> genmap_lmers_starting_positions, int no_of_sequences, int sequence_length, StringSet<DnaString> sequences) { //added &
    cout << "running genmap2 ..." << endl;
    cout << "size of genmap_lmers_starting_positions is " << genmap_lmers_starting_positions.size() << endl;
    //convert every starting position into a bucket so this bucket can be used to create the initial weight matrix initWh
    //for loop that iterates through genmap_lmers_starting_positions
    // for each pair in the vector of pairs genmap_lmers_starting_positions, put the pair in the function convertPairToMap and call for (bucket : pair_as_map)
    pair<DnaString, int> best_genmap_conseq;
    int best_genmap_score = -1;
    for(int i = 0; i < genmap_lmers_starting_positions.size(); i++) {
        map<int, vector<pair<int, int>>> pair_as_map = convertPairToMap(genmap_lmers_starting_positions.at(i));
        vector<pair<DnaString, int>> genmap_bucket_conseqs;
        for (pair<int, vector<pair<int, int>>> bucket: pair_as_map) {
            vector<vector<float>> Wh;
            vector<vector<float>> posM = vector<vector<float>>(length(sequences), (vector<float>(length(sequences[0]) - motif_length + 1, 0)));
            initWh(motif_length, bucket, Wh); //initialize weight matrix Wh for prob of a base in the motif
            for (int refine_iter = 0; refine_iter < 10; refine_iter++) //Refine weight matrix W and posM until convergence
                refine(motif_length, Wh, posM);
            get_consensus_seq(motif_length, d, posM, genmap_bucket_conseqs); //CONSENSUS SEQUENCE
            //cout << "getting here" << endl;
            pair<DnaString, int> curr_genmap_conseq = genmap_bucket_conseqs.at(0);
            int curr_genmap_score = curr_genmap_conseq.second;
            if (best_genmap_score == -1 || curr_genmap_score <= best_genmap_score) {
                best_genmap_score = curr_genmap_score;
                best_genmap_conseq.first = curr_genmap_conseq.first;
                best_genmap_conseq.second = curr_genmap_conseq.second;
            }
        }
        //cout << "getting there" << endl;
    }
    double lmer_performance_coefficient = 0;
    int exact_matches = 0;
    //double genmap_lmer_performance_coefficient = calculatePerformanceCoefficient(best_genmap_conseq.first, lmer_performance_coefficient, exact_matches);
    cout << "getting here too" << endl;
    cout << "searched consensus sequence: [" << best_genmap_conseq.first << "] with a score of " << best_genmap_conseq.second << endl;
    //cout << "The performance coefficient is " << roundWithPrecision(genmap_lmer_performance_coefficient, 2) << "." << endl;
    return best_genmap_conseq.first;
};

//gm

/**
 * Makes use of the seqan::find function to find the computed pattern in the sequences.
 * We want to print out the found matches and their positions.
 */
void printPositions(int num_of_seqs, DnaString pattern) {
    auto delegate = [](auto & iter, DnaString const & needle, uint8_t errors) {
        for (auto occ : getOccurrences(iter))
            addToMap<pair<int,DnaString>>(foundMatches, seqidx, pair<int,DnaString>(occ,substr(sequences[seqidx], occ, l)));
    };

    for(seqidx=0; seqidx<num_of_seqs; seqidx++) {
        Index<DnaString, BidirectionalIndex<FMIndex<> > > index(sequences[seqidx]);
        switch (d) {
            case 1: find<0,1>(delegate, index, pattern, HammingDistance()); break;
            case 2: find<0,2>(delegate, index, pattern, HammingDistance()); break;
            case 3: find<0,3>(delegate, index, pattern, HammingDistance()); break;
            case 4: find<0,4>(delegate, index, pattern, HammingDistance()); break;
        }
    }

    for(auto p : foundMatches) { // print out the foundMatches
        cout << p.first << ": { ";
        for(pair<int,DnaString> pos : p.second)
            cout << "[" << pos.first << " => " << pos.second << "] ";
        cout << "}" << endl;
    }
}



/**
 * If a planted motif csv is given, check the found matches against the
 * positions in the planted motif file. Compute then a performance coefficient and print it.ca
 */
void printPerformanceCoefficient(DnaString pattern) {
    int z=0;
    int n=0;

    if(DnaString(pm[0][2]) == pattern) {
        z=1; n=1; exactMatches++;
    } else {

        for(int i=0; i<length(sequences); i++) {
            int planted_position = stoi(pm[i+1][1]);
            if(foundMatches.count(i)>0) {
                int lowerBound = planted_position-(l-1);
                int upperBound = planted_position+(l-1);
                int u = l;
                int current_fm;
                for(pair<int,DnaString> foundMatch: foundMatches[i]) {
                    current_fm = foundMatch.first;
                    if(current_fm >= lowerBound && current_fm <= upperBound)
                        break;
                }
                int tail = (current_fm >= lowerBound && current_fm <= upperBound) ? abs(planted_position-current_fm) : l;
                u = min(u,tail);
                z += (l - u);
                n += (l + u);
            } else {
                n+=l;
            }
        }
    }
    avrg_performance_coefficient += (double)z/(double)n;
    cout << "The performance coefficient is " << rwp((double)z/(double)n,2) << "." << endl;
}



int main(int argc, char const ** argv) {
    program_start = chrono::steady_clock::now();
    noOfThreads = thread::hardware_concurrency();
    if (argc < 3) {
        cerr << "Not enough arguments, expected are at least 2: genome-file (in fasta) and parameters file (in csv). "
                "Optional parameters are: planted motiff file (in csv) and amount of datasets (integer).\n";
        return -1;
    }
    string fastaFilename(argv[1]);



    //gm
    size_t lastSlash = fastaFilename.find_last_of("/");
    string file_path = fastaFilename.substr(0, lastSlash);
    string filenameOnly = fastaFilename.substr(lastSlash+1);
    size_t ftpos0 = filenameOnly.find(".fasta");
    string file_without_suffix = filenameOnly.substr(0, ftpos0);
    //gm





    size_t ftpos = fastaFilename.find(".fasta");
    if(ftpos == -1) {cerr << "invalid fasta file. Filename has to end with .fasta\n"; return -1;}
    string datasetName = fastaFilename.substr(0,ftpos);
    string parametersFile(argv[2]), plantedMotifFilename, pmName;
    if(argc > 3) {
        plantedMotifFilename = string(argv[3]);
        ftpos = plantedMotifFilename.find(".csv");
        if(ftpos == -1) {cerr << "invalid planted motif file. Filename has to end with .csv\n"; return -1;}
        pmName = plantedMotifFilename.substr(0,ftpos);
    }


    ifstream file2(parametersFile);
    if(file2.good()) parameters = csvParseIntoVector(file2,'\t',1);
    else { cerr << "Parameters file not found or not readable\n"; return -1; }

    cout << endl;
    l=stoi(parameters[0][0]);
    d=stoi(parameters[0][1]);/*allowed mutations*/
    k=stoi(parameters[0][2]);/*the number of positions to be projected*/
    s=stoi(parameters[0][3]);//choose a value for s (bucket threshold)
    m=stoi(parameters[0][4]);//calculate the optimal number of trials m;
    datasets = 1;

    if(argc == 5) datasets = stoi(argv[4]);
    string no = "";

    datasets_start = chrono::steady_clock::now();
    step = datasets_start;
    for(int i = 1; i<=datasets; i++) {
        tmp_round = i;
        chrono::steady_clock::time_point dataset_start = chrono::steady_clock::now();
        if(argc == 5) no = "_"+to_string(i);

        //use the fai file to get the content of the fasta file
        FaiIndex faiIndex;
        if (!build(faiIndex, (datasetName+no+".fasta").c_str())) {
            cerr << "ERROR: Could not build FAI index for file "<< (datasetName+no+".fasta").c_str() << "\n";
            return -1;
        }

        if (!save(faiIndex, toCString((datasetName+no+".fasta.fai").c_str()))) {
            cerr << "ERROR: Could not write the index to file!\n";
            return -1;
        }

        cout << "Index file " << (datasetName+no+".fasta.fai").c_str() << " was successfully created.\n";

        if (!open(faiIndex, (datasetName+no+".fasta").c_str()))
            cout << "ERROR: Could not load FAI index " << (datasetName+no+".fasta.fai").c_str() << "\n";

        chrono::steady_clock::time_point start = chrono::steady_clock::now();

        //gm
        /*
        * Compute the mappbility/get genmap frequency vector
        */
        vector<uint8_t> genmap_frequency_vector = getGenMapFrequencyVector(datasetName + no, file_without_suffix + no, l, d);
        //gm
        if(argc>3) {
            ifstream file(pmName+no+".csv");
            if(file.good()) pm = csvParseIntoVector(file,'\t',1);
            else { cerr << "Planted motif file not found or not readable\n"; return -1; }
        }
        DnaString pattern;
        unsigned num_of_seqs = numSeqs(faiIndex);
        for(unsigned idx = 0; idx < num_of_seqs; idx++){
            DnaString seq_in_file;
            readSequence(seq_in_file, faiIndex, idx);
            appendValue(sequences, seq_in_file); //save each sequence in the stringSet sequences
            //for genmap create a fasta-file for every sequence in the dataset

                string mkdirr = "mkdir -p ./genmap_fasta_files/genmap_" + file_without_suffix + no;
                const int dir_error = system(mkdirr.c_str());
                if (-1 == dir_error)
                {
                    printf("Error creating subdirectory!n");
                    exit(1);
                }
                string fasta_file_line_name = "./genmap_fasta_files/genmap_" + file_without_suffix + no + "/" + file_without_suffix + no + "_" + to_string(idx+1) + ".fasta";
                cout << fasta_file_line_name << endl;
                ofstream file(fasta_file_line_name);
                file << ">seq0" << "\n" << seq_in_file;

        }
        //for each dataset (eg folder with files that contain only one sequence) use genmap on a folder with --exclude-pseudo in order to count the number of files that contain an lmer with up to d mismatches
        string folder_with_fasta_file_lines = "./genmap_fasta_files/genmap_" + file_without_suffix + no;
        cout << "file_without_suffix is " << file_without_suffix << endl;
        vector<pair<int,int>> lmers_contained_in_many_files;
        int min_no_of_files = num_of_seqs*0.5;

            for(unsigned idx = 0; idx < num_of_seqs; idx++){
                string folder_with_fasta_file_lines_filename = "genmap_" + file_without_suffix + no + "_" + to_string(idx+1);
                vector<uint8_t> frequency_vector_freq8 = getGenMapFrequencyVectorOPS(folder_with_fasta_file_lines, folder_with_fasta_file_lines_filename, l, d);
                vector<int> frequency_vector_int = covertFreqVecToIntVec(frequency_vector_freq8);
                //saveLmerIfInMinNoOfFiles(frequency_vector_freq8, lmers_contained_in_many_files, idx, min_no_of_files);
                for (int j = 0; j < length(frequency_vector_freq8); j++) {
                    if(frequency_vector_freq8[j] >= min_no_of_files) {
                        lmers_contained_in_many_files.push_back(pair<int,int>(idx,j));
                    }
                }
            }



        pattern = runGenMap2(l, lmers_contained_in_many_files, length(sequences), length(sequences[0]), sequences);

        //gm


        if(length(pattern) == 0) return -4; // No sequences in any bucket
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        times.push_back(chrono::duration_cast<chrono::milliseconds>(end-start).count()/m);

        if(d<5) printPositions(num_of_seqs, pattern);
        if(argc>3) printPerformanceCoefficient(pattern);

        done_datasets++; // One dataset is processed
        clear(sequences); foundMatches.clear();	seqidx=0; // reset
        chrono::steady_clock::time_point dataset_end = chrono::steady_clock::now();
        cout << "Processing this dataset took " << chrono::duration_cast<chrono::milliseconds>(dataset_end-dataset_start).count()/1000 << " seconds." << endl<<endl;
    }


    double average_time;
    for(double time : times) average_time+=time;
    average_time = average_time/length(times);
    cout << "The average time for one trial was " << average_time/1000 << " seconds." << endl;

    cout << "For " << exactMatches << " out of " << datasets << " dataset(s) the correct planted motifs were found. The average performance coefficient is " << rwp(avrg_performance_coefficient/(double) datasets,2) << endl;

    chrono::steady_clock::time_point program_end = chrono::steady_clock::now();
    vector<int> time_passed = secondsToHours(chrono::duration_cast<chrono::seconds>(program_end-program_start).count());
    cout << "Execution of program finished after (hh:mm:ss) " << timeVectorToString(time_passed) << "." << endl;
    return 0; // execution successful
}
