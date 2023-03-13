#include <iostream>
#include <chrono>
#include <cstdlib>

#include "br_index.hpp"
#include "br_index_nplcp.hpp"
#include "utils.hpp"

using namespace bri;
using namespace std;

ulint kappa = 1; // MEM threshold
long allowed = 0;
bool nplcp = false;
bool relax = false;
string alphabet = "";
string output_file = "";

void help()
{
	cout << "bri-mem: locate all MEMs between patterns and text" << endl;

	cout << "Usage: bri-mem [options] <index> <patterns> [<output>]" << endl;
    cout << "   -nplcp       use the version without PLCP " << endl;
    cout << "   -relax    use relaxed definition for MEMs " << endl;
	cout << "   -k      MEM threshold" << endl;
	cout << "   -a      alphabet string with last symbol regarded as separator, default ACGT#" << endl;	
	cout << "   -o      output file where lines i,x,d are outputed for MEMs P[i..i+d-1]=T[x..x+d-1] " << endl;
	cout << "   -af     file containing alphabet " << endl;
	cout << "   <text>      text index file (with extension .bri)" << endl;
	cout << "   <patterns>  index file of concatenation of patterns (with extension .bri)" << endl;
	cout << "               concatenation format: #AGGATG#AGATGT#, where # is separator symbol" << endl;		
	exit(0);
}

bool parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

    if (s.compare("-nplcp") == 0)
    {

        nplcp = true;

    }
    if (s.compare("-relax") == 0)
    {

        relax = true;
        cout << "Relaxed MEM reporting not yet supported." << endl;

    }
	else if (s.compare("-k") == 0) 
    {

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -k option." << endl;
			help();
		}

		kappa = atoi(argv[ptr]);
		ptr++;
    
    }
    else if (s.compare("-a") == 0)
    {

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -a option." << endl;
			help();
		}
        alphabet = string(argv[ptr]);
        ptr++;
    }
    else if (s.compare("-o") == 0)
    {

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -o option." << endl;
			help();
		}
        output_file = string(argv[ptr]);  
        ptr++;
    }
    else if (s.compare("-af") == 0)
    {

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -af option." << endl;
			help();
		}
        string alphabet_file="";
        alphabet_file = string(argv[ptr]);
        ifstream as(alphabet_file);
        getline(as,alphabet);
        as.close();    
    }
    else
    {
        return 0;
	}
	return 1;

}

template<class T, class TS>
void reportMEMs(T pidx, T idx, TS psample, TS sample, ulint d, ofstream& output)
{
   // a bit naive implementation of the cross product 
   std::vector<ulint>** a= new  std::vector<ulint>*[alphabet.size()];
   std::vector<ulint>** b= new  std::vector<ulint>*[alphabet.size()];
   TS left,right;
   for (ulint i=0; i<alphabet.size(); i++) {
      a[i] = new std::vector<ulint>[alphabet.size()];
      b[i] = new std::vector<ulint>[alphabet.size()];      
      for (ulint j=0; j<alphabet.size(); j++) {
         right = pidx.right_extension(alphabet[j],psample);
         if (!right.is_invalid()) {
            left = pidx.left_extension(alphabet[i],right);
            if (!left.is_invalid())
               a[i][j] = pidx.locate_sample(left);
         }   
         right = idx.right_extension(alphabet[j],sample);
         if (!right.is_invalid()) {
            left = idx.left_extension(alphabet[i],right);
            if (!left.is_invalid())
               b[i][j] = idx.locate_sample(left);
         }
      }
   }
   for (ulint i=0; i<alphabet.size(); i++) 
      for (ulint j=0; j<alphabet.size(); j++)
         if (a[i][j].size()>0) {
            for (ulint ii=0; ii<alphabet.size(); ii++)
               for (ulint jj=0; jj<alphabet.size(); jj++)
                  if (b[ii][jj].size()>0 and i!=ii and j!=jj) {
                     for (ulint iii=0; iii<a[i][j].size(); iii++) 
                        for (ulint jjj=0; jjj<b[ii][jj].size(); jjj++)
                            output << a[i][j][iii]+1 << "," << b[ii][jj][jjj]+1 << "," << d << endl;          
                  }
         }
   for (ulint i=0; i<alphabet.size(); i++) {
      delete[] a[i];
      delete[] b[i];
   }
   delete[] a;
   delete[] b;    
}

template<class T,class TS>
void find_mems(ifstream& in, ifstream& pin)
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    using std::chrono::microseconds;


    auto t1 = high_resolution_clock::now();

    T idx;

    idx.load(in);
    
    T pidx;
    pidx.load(pin);

    cout << "Searching MEMs " << endl;
    
    TS sample(idx.get_initial_sample(true));
    TS psample(pidx.get_initial_sample(true));
    pair <TS,TS> node;
    std::stack<pair <TS,TS>> S; // interval pairs
    std::stack<ulint> dS; // string depths
    node.first = sample;
    node.second = psample;
    // node is now suffix tree root
    S.push(node);
    ulint d = 0;
    dS.push(d); 
    ulint x;
    bool MEM;
    ulint maxMEM = 0;
    
    if (alphabet.size()==0) {
       alphabet = "ACGT#";
    }
    
    ofstream output;
    if (output_file.size()!=0) {     
       output.open(output_file);
       if (!output) {
          cout << "Could not open output file " << output_file << endl;
       }
    }
    
    while (!S.empty()) {
       node = S.top();
       S.pop();
       d = dS.top();
       dS.pop();
       sample = node.first;
       psample = node.second;
       if (sample.is_invalid() or 
          psample.is_invalid()) { 
          continue; // not a valid range
       }
       if ((!idx.is_right_maximal(sample) and !pidx.is_right_maximal(psample)) and 
           idx.bwt_at(sample.rangeR.first,true)==pidx.bwt_at(psample.rangeR.first,true) ) {
          continue; // implicit node reached
       }

       MEM = 1;
       // Taking Weiner links from current node to visit nodes at string depth++
       // TODO: Push the largest interval first to limit the size of the stack
       // TODO: use as alphabet the distinct symbols appearing in the range

       // last alphabet symbol regarded as separator
       for (ulint j=0;j<alphabet.size()-1; j++) { 
          node.first = idx.left_extension(alphabet[j],sample); 
          node.second = pidx.left_extension(alphabet[j],psample);   
          S.push(node);       
          dS.push(d+1);
          // range not splitting, no MEM reported
          if (sample.size()+psample.size()==node.first.size()+node.second.size()) 
             MEM = 0;
       }
       if (d > maxMEM)
          maxMEM = d;
       if (MEM and d>=kappa and output.is_open())
             reportMEMs<T,TS>(pidx,idx,psample,sample,d,output);    
    }
    cout << "Maximum MEM is of length " << maxMEM << endl;

    output.close();

    auto t2 = high_resolution_clock::now();
    ulint tot_time = duration_cast<microseconds>(t2-t1).count(); 
    cout << "Total time     : " << tot_time << " microseconds" << endl;
}



int main(int argc, char** argv)
{
    if (argc < 3) help();

    int ptr = 1;

    while (ptr < argc - 2) 
       parse_args(argv, argc, ptr);

    string idx_file(argv[ptr]);
    string patt_file(argv[ptr+1]);
    
    cout << "Loading br-indexes" << endl;

    ifstream in(idx_file);
    ifstream pin(patt_file);
    
    if (nplcp) {
       find_mems<br_index_nplcp<>,br_sample_nplcp >(in, pin);
    }
    else {
       find_mems<br_index<>,br_sample >(in, pin);
    }

    in.close();
    pin.close();

}
