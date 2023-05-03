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
bool asymmetric = false;
string alphabet = "";
string output_file = "";

void help()
{
	cout << "bri-mem: locate all MEMs between queries and text" << endl;

	cout << "Usage: bri-mem [options] <index> <queries>" << endl;
        cout << "   --nplcp       use the version without PLCP " << endl;
        cout << "   --asymmetric  use asymmetric definition for MEMs " << endl;
	cout << "   -k      MEM threshold" << endl;
	cout << "   -a      alphabet string with last symbol regarded as separator, default ACGT#" << endl;	
	cout << "   -o      output file where lines x,i,d are outputed for MEMs Q[i..i+d-1]=T[x..x+d-1]; " << endl;
        cout << "           in symmetric mode (default) the matches are locally maximal, i.e., Q[i-1]!=T[x-1] and Q[i+d]!=T[x+d] " << endl;
	cout << "           in asymmetric mode the matches are globally maximal, i.e., Q[i..i+d] or Q[i-1..i+d-1] do not " << endl;
	cout << "           occur in T; only one occurrence T[x..x+d-1] is reported in asymmetric mode" << endl;		
	cout << "   -f     file containing alphabet " << endl;
	cout << "   <text>      text index file (with extension .bri)" << endl;
	cout << "   <queries>  index file of concatenation of queries (with extension .bri)" << endl;
	cout << "               concatenation format: #AGGATG#AGATGT#, where # is separator symbol" << endl;		
	exit(0);
}

bool parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

    if (s.compare("--nplcp") == 0)
    {

        nplcp = true;

    }
    if (s.compare("--asymmetric") == 0)
    {

        asymmetric = true;

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
    else if (s.compare("-f") == 0)
    {

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -f option." << endl;
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
void reportMEMs(T qidx, T idx, TS qsample, TS sample, ulint d, ofstream& output)
{
   // a bit naive implementation of the cross product:
   // factor alphabet.size()^2 slower than an optimal implementation
    
   std::vector<ulint>** a= new  std::vector<ulint>*[alphabet.size()];
   std::vector<ulint>** b= new  std::vector<ulint>*[alphabet.size()];
   TS** Sa= new  TS*[alphabet.size()];
   TS** Sb= new  TS*[alphabet.size()];
   bool** Ia= new  bool*[alphabet.size()];
   bool** Ib= new  bool*[alphabet.size()];
 
   TS left,right;
   for (ulint i=0; i<alphabet.size(); i++) {
      a[i] = new std::vector<ulint>[alphabet.size()];
      b[i] = new std::vector<ulint>[alphabet.size()];      
      Sa[i] = new TS[alphabet.size()];
      Sb[i] = new TS[alphabet.size()];
      Ia[i] = new bool[alphabet.size()];
      Ib[i] = new bool[alphabet.size()];      
      for (ulint j=0; j<alphabet.size(); j++) {
         Ia[i][j] = false;
         right = qidx.right_extension(alphabet[j],qsample);
         if (!right.is_invalid()) {
            left = qidx.left_extension(alphabet[i],right);         
            if (!left.is_invalid()) {
               Sa[i][j] = left;
               Ia[i][j] = true;
            }
         }   
         Ib[i][j] = false;
         right = idx.right_extension(alphabet[j],sample);
         if (!right.is_invalid()) {
            left = idx.left_extension(alphabet[i],right);
            if (!left.is_invalid()) {
               Sb[i][j] = left;
               Ib[i][j] = true;
            }
         }
      }
   }
   for (ulint i=0; i<alphabet.size(); i++) 
      for (ulint j=0; j<alphabet.size(); j++)
         if (Ia[i][j]) {
            for (ulint ii=0; ii<alphabet.size(); ii++)
               for (ulint jj=0; jj<alphabet.size(); jj++)
                  if (Ib[ii][jj] and i!=ii and j!=jj) {
                     if (a[i][j].size()==0)  // locate charged on the output
                        a[i][j] = qidx.locate_sample(Sa[i][j]);
                     if (b[ii][jj].size()==0)  // locate charged on the output
                        b[ii][jj] = idx.locate_sample(Sb[ii][jj]);
                     for (ulint iii=0; iii<a[i][j].size(); iii++) 
                        for (ulint jjj=0; jjj<b[ii][jj].size(); jjj++)
                            output <<  b[ii][jj][jjj]+1 << "," << a[i][j][iii]+1 << "," << d << endl;          
                  }
         }
   for (ulint i=0; i<alphabet.size(); i++) {
      delete[] a[i];
      delete[] b[i];
      delete[] Sa[i];
      delete[] Sb[i];
      delete[] Ia[i];
      delete[] Ib[i];      
   }
   delete[] a;
   delete[] b;    
   delete[] Sa;
   delete[] Sb;    
   delete[] Ia;
   delete[] Ib;    
}

template<class T, class TS>
void reportAMEMs(T qidx, T idx, TS qsample, TS sample, ulint d, ofstream& output)
{
   std::vector<ulint>** a= new  std::vector<ulint>*[alphabet.size()];
   bool** b= new  bool*[alphabet.size()];
   TS** Sa= new  TS*[alphabet.size()];
   bool** Ia= new  bool*[alphabet.size()];

   TS left,right;
   for (ulint i=0; i<alphabet.size(); i++) {
      a[i] = new std::vector<ulint>[alphabet.size()];
      b[i] = new bool[alphabet.size()];
      Sa[i] = new TS[alphabet.size()];
      Ia[i] = new bool[alphabet.size()];
      for (ulint j=0; j<alphabet.size(); j++) {
         Ia[i][j] = false;
         right = qidx.right_extension(alphabet[j],qsample);
         if (!right.is_invalid()) {
            left = qidx.left_extension(alphabet[i],right);
            if (!left.is_invalid()) {
               Sa[i][j] = left;
               Ia[i][j] = true;
            }
         }   
         b[i][j] = true; // maximal in text?         
         right = idx.right_extension(alphabet[j],sample);
         if (!right.is_invalid())
            b[i][j] = false; // not right-maximal  
         left = idx.left_extension(alphabet[i],sample);
         if (!left.is_invalid())
            b[i][j] = false; // not left-maximal
      }
   }
   for (ulint i=0; i<alphabet.size(); i++) 
      for (ulint j=0; j<alphabet.size(); j++)
         if (Ia[i][j] and b[i][j]) { 
            a[i][j] = qidx.locate_sample(Sa[i][j]);
            for (ulint iii=0; iii<a[i][j].size(); iii++) 
               output << a[i][j][iii]+1 << ","  << sample.j -sample.d << "," << d << endl;          
         }
   for (ulint i=0; i<alphabet.size(); i++) {
      delete[] a[i];
      delete[] b[i];
      delete[] Sa[i];
      delete[] Ia[i];
      
   }
   delete[] a;
   delete[] b;    
   delete[] Sa;
   delete[] Ia;
}

template<class T,class TS>
void find_mems(ifstream& in, ifstream& qin)
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    using std::chrono::microseconds;


    auto t1 = high_resolution_clock::now();

    T idx;

    idx.load(in);
    
    T qidx;
    qidx.load(qin);

    cout << "Searching MEMs " << endl;
    
    TS sample(idx.get_initial_sample(true));
    TS qsample(qidx.get_initial_sample(true));
    pair <TS,TS> node;
    std::stack<pair <TS,TS>> S; // interval pairs
    std::stack<ulint> dS; // string depths
    node.first = sample;
    node.second = qsample;
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
       qsample = node.second;
       if (sample.is_invalid() or 
          qsample.is_invalid()) { 
          continue; // not a valid range
       }
       if ((!idx.is_right_maximal(sample) and !qidx.is_right_maximal(qsample)) and 
           idx.bwt_at(sample.rangeR.first,true)==qidx.bwt_at(qsample.rangeR.first,true) ) {
          continue; // implicit node reached
       }

       MEM = 1;
       // Taking Weiner links from current node to visit nodes at string depth++
       // TODO: Push the largest interval first to limit the size of the stack
       // TODO: use as alphabet the distinct symbols appearing in the range

       // last alphabet symbol regarded as separator
       for (ulint j=0;j<alphabet.size()-1; j++) { 
          node.first = idx.left_extension(alphabet[j],sample); 
          node.second = qidx.left_extension(alphabet[j],qsample);   
          S.push(node);       
          dS.push(d+1);
          // range not splitting, no MEM reported
          if (sample.size()+qsample.size()==node.first.size()+node.second.size()) 
             MEM = 0;
       }
       if (d > maxMEM)
          maxMEM = d;
       if (MEM and d>=kappa and output.is_open())
          if (!asymmetric)
             reportMEMs<T,TS>(qidx,idx,qsample,sample,d,output);
          else 
             reportAMEMs<T,TS>(qidx,idx,qsample,sample,d,output);   
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
    string query_file(argv[ptr+1]);
    
    cout << "Loading br-indexes" << endl;

    ifstream in(idx_file);
    ifstream qin(query_file);
    
    if (nplcp) {
       find_mems<br_index_nplcp<>,br_sample_nplcp >(in, qin);
    }
    else {
       find_mems<br_index<>,br_sample >(in, qin);
    }

    in.close();
    qin.close();

}
