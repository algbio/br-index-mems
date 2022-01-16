#include <iostream>
#include <chrono>
#include <cstdlib>

#include "br_index.hpp"
#include "utils.hpp"

using namespace bri;
using namespace std;

string check = string();
long allowed = 0;

void help()
{
	cout << "bri-locate: locate all occurrences of the input patterns" << endl;
    cout << "             allowing some mismatched characters."        << endl << endl;

	cout << "Usage: bri-locate [options] <index> <patterns>" << endl;
    cout << "   -m <number>  max number of mismatched characters allowed (supported: 0,1,2 (0 by default))" << endl;
	cout << "   -c <text>    check correctness of each pattern occurrence on this text file (must be the same indexed)" << endl;
	cout << "   <index>      index file (with extension .bri)" << endl;
	cout << "   <patterns>   file in pizza&chili format containing the patterns." << endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if(s.compare("-c")==0){

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -c option." << endl;
			help();
		}

		check = string(argv[ptr]);
		ptr++;
    
    }else if(s.compare("-m")==0){

        if(ptr>=argc-1){
            cout << "Error: missing parameter after -m option." << endl;
            help();
        }

        char* e;
        allowed = strtol(argv[ptr],&e,10);

        if(*e != '\0' || allowed < 0 || allowed > 2){
            cout << "Error: invalid value not 0, 1, 2 after -m option." << endl;
            help();
        }

        ptr++;

	}else{

		cout << "Error: unknown option " << s << endl;
		help();

	}

}


void locate_all(ifstream& in, string patterns)
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    string text;
    bool c = false;

    if (check.compare(string()) != 0)
    {
        c = true;

        ifstream ifs1(check);
        stringstream ss;
        ss << ifs1.rdbuf();
        text = ss.str();
    }

    auto t1 = high_resolution_clock::now();

    br_index<> idx;

    idx.load(in);

    auto t2 = high_resolution_clock::now();

    cout << "searching patterns with mismatches at most " << allowed << " ... " << endl;
    ifstream ifs(patterns);

    //read header of the pizza&chilli input file
	//header example:
	//# number=7 length=10 file=genome.fasta forbidden=\n\t
    string header;
    getline(ifs,header);

    ulint n = get_number_of_patterns(header);
    ulint m = get_patterns_length(header);

    ulint last_perc = 0;

    ulint occ0 = 0, occ1 = 0, occ2 = 0;
    ulint occ_tot = 0;

    // extract patterns from file and search them in the index
    for (ulint i = 0; i < n; ++i)
    {
        ulint perc = 100 * i / n;
        if (perc > last_perc)
        {
            cout << perc << "% done ..." << endl;
            last_perc = perc;
        }

        string p = string();

        for (ulint j = 0; j < m; ++j)
        {
            char c;
            ifs.get(c);
            p += c;
        }

        // exact match
        vector<ulint> occs = idx.locate(p);
        occ0 += occs.size();
        occ_tot += occs.size();

        if (c) // check occurrences
        {
            cout << "number of exact occs             : " << occs.size() << endl;
            for (auto o : occs)
            {
                if (text.substr(o,p.size()).compare(p) != 0) 
                {
                    cout << "Error: wrong occurrence:  " << o << endl;
                    cout << "       original pattern:  " << p << endl;
                    cout << "       wrong    pattern: "  << text.substr(o,p.size()) << endl;
                }
            }
        }

        if (allowed < 1) continue;

        // approximate match with 1 miss
        occs = idx.locate1(p);
        occ1 += occs.size();
        occ_tot += occs.size();

        if (c) // check occurrences
        {
            cout << "number of occs with 1 mismatch   : " << occs.size() << endl;
            for (auto o : occs)
            {
                int mismatches = 0;
                for (size_t i = 0; i < p.size(); ++i)
                {
                    if (text[o+i] != p[i]) mismatches++;
                }
                if (mismatches != 1) 
                {
                    cout << "Error: wrong occurrence:  " << o << endl;
                    cout << "       original pattern:  " << p << endl;
                    cout << "       wrong    pattern: " << text.substr(o,p.size()) << endl;
                }
            }
        }

        if (allowed < 2) continue;

        // approximate match with 2 miss
        occs = idx.locate2(p);
        occ2 += occs.size();
        occ_tot += occs.size();

        if (c) // check occurrences
        {
            cout << "number of occs with 2 mismatches : " << occs.size() << endl;
            for (auto o : occs)
            {
                int mismatches = 0;
                for (size_t i = 0; i < p.size(); ++i)
                {
                    if (text[o+i] != p[i]) mismatches++;
                }
                if (mismatches != 2) 
                {
                    cout << "Error: wrong occurrence:  " << o << endl;
                    cout << "       original pattern:  " << p << endl;
                    cout << "       wrong    pattern:  " << text.substr(o,p.size()) << endl;
                }
            }
        }
    }

    double occ_avg = (double)occ_tot / n;
    
    cout << endl << occ_avg << " average occurrences per pattern" << endl;

    ifs.close();

    auto t3 = high_resolution_clock::now();

    ulint load = duration_cast<milliseconds>(t2-t1).count();
    cout << "Load time  : " << load << " milliseconds" << endl;

    ulint search = duration_cast<milliseconds>(t3-t2).count();
    cout << "Number of patterns n = " << n << endl;
	cout << "Pattern length     m = " << m << endl;
    cout << "Occurrences with 0 mismatch   occ0 = " << occ0 << endl;
    if (allowed >= 1)
    {
        cout << "Occurrences with 1 mismatch   occ1 = " << occ1 << endl;
        if (allowed >= 2)
        {
            cout << "Occurrences with 2 mismatches occ2 = " << occ2 << endl;
        }
    }
	cout << "Total number of occurrences   occt = " << occ_tot << endl << endl;

    cout << "Total time : " << search << " milliseconds" << endl;
	cout << "Search time: " << (double)search/n << " milliseconds/pattern (total: " << n << " patterns)" << endl;
	cout << "Search time: " << (double)search/occ_tot << " milliseconds/occurrence (total: " << occ_tot << " occurrences)" << endl;
}



int main(int argc, char** argv)
{
    if (argc < 3) help();

    int ptr = 1;

    while (ptr < argc - 2) parse_args(argv, argc, ptr);

    string idx_file(argv[ptr]);
    string patt_file(argv[ptr+1]);

    ifstream in(idx_file);

    cout << "Loading br-index" << endl;
    locate_all(in, patt_file);

    in.close();

}