#include <iostream>
#include <chrono>
#include <cstdlib>

#include <sdsl/bit_vectors.hpp>
#include "br_index.hpp"
#include "br_index_nplcp.hpp"
#include "utils.hpp"

using namespace bri;
using namespace std;
using namespace sdsl;

// Global variables assigned with optional parameters
ulint kappa = 1; // MEM threshold
bool nplcp = false;
bool asymmetric = false;
bool filter = false;
string alphabet = "";
string output_file = "";

void help()
{
	cout << "bri-mem-efg: locate all MEMs between queries and elastic founder graph" << endl;

	cout << "Usage: bri-mem-efg [options] <text> <queries> <efg>" << endl;
    cout << "   --nplcp       use the version without PLCP " << endl;
    cout << "   --asymmetric  use asymmetric definition for MEMs " << endl;
	cout << "   --filter      filter out MEMs that do not occur in <text> " << endl;
	cout << "   -k      MEM threshold" << endl;
	cout << "   -a      alphabet string with last symbol regarded as separator, default ACGT#" << endl;	
	cout << "   -o      output file where lines x,i,d are outputed for MEMs Q[i..i+d-1]=T[x..x+d-1]; " << endl;
    cout << "           in symmetric mode (default) the matches are locally maximal, i.e., Q[i-1]!=T[x-1] and Q[i+d]!=T[x+d] " << endl;
	cout << "           in asymmetric mode the matches are globally maximal, i.e., Q[i..i+d] or Q[i-1..i+d-1] do not " << endl;
	cout << "           occur in T; only one occurrence T[x..x+d-1] is reported in asymmetric mode" << endl;		
	cout << "           output is formatted as fasta file with MEMs listed according to the input efg fasta" << endl;		
	cout << "   -f     file containing alphabet " << endl;
	cout << "   <text>  index file of the concatenation of MSA rows (with extension .bri)" << endl;
	cout << "   <queries>  index file of concatenation of queries (with extension .bri)" << endl;
	cout << "               concatenation format: #AGGATG#AGATGT#, where # is separator symbol" << endl;		
	cout << "   <efg>      fasta file with labels of paths of length 1, 2, and 3 concatenated using gsa2concats < .gfa" << endl;
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
    else if (s.compare("--asymmetric") == 0)
    {

        asymmetric = true;

    }
    else if (s.compare("--filter") == 0)
    {

        filter = true;

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
void reportMEMs(T idx, T qidx, TS sample, TS qsample, ulint d, ofstream& output, bit_vector Bsuffix = bit_vector(), bit_vector Bprefix = bit_vector())
{ 
   // a bit naive implementation of the cross product 
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
                  // The following is outputting some non-maximal matches, as we don't check if a MEM has unique left or right context in the graph
                  if (Ib[ii][jj] and (i!=ii or (i==alphabet.size()-1 and ii==alphabet.size()-1))  and (j!=jj or (j==alphabet.size()-1 and jj==alphabet.size()-1))) {
                     if (a[i][j].size()==0)  // locate charged on the output
                        a[i][j] = qidx.locate_sample(Sa[i][j]);
                     if (b[ii][jj].size()==0)  // locate charged on the output
                        b[ii][jj] = idx.locate_sample(Sb[ii][jj]);                  
                     for (ulint iii=0; iii<a[i][j].size(); iii++) 
                        for (ulint jjj=0; jjj<b[ii][jj].size(); jjj++)
                           if (Bsuffix.size()==0 or Bprefix.size()==0 or (Bsuffix[b[ii][jj][jjj]+1] and Bprefix[b[ii][jj][jjj]+d]))
                              output << b[ii][jj][jjj]+1 << "," << a[i][j][iii]+1 << "," << d << endl;          
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

template<class T, class TS>
ulint explore_mems(T tidx, T qidx, T fidx, ofstream& output, bool f = false, bit_vector Bsuffix = bit_vector(), bit_vector Bprefix = bit_vector())
{
    TS sample(tidx.get_initial_sample(true));
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
    bool MEM;
    ulint maxMEM = 0;    
    
    TS fsample;
    if (f) // using fidx as filter
       fsample = fidx.get_initial_sample(true);
         
    std::stack<TS> fS; // filter text index range
    if (f)
       fS.push(fsample);
     
    while (!S.empty()) {
       node = S.top();
       S.pop();
       d = dS.top();
       dS.pop();
       sample = node.first;
       qsample = node.second;
       if (f) {
          fsample = fS.top();
          fS.pop();
       }
       
       if (sample.is_invalid() or 
          qsample.is_invalid() or
          (!asymmetric and f and fsample.is_invalid())) { 
          continue; // not a valid range
       }
       if ((!tidx.is_right_maximal(sample) and !qidx.is_right_maximal(qsample)) and 
           tidx.bwt_at(sample.rangeR.first,true)==qidx.bwt_at(qsample.rangeR.first,true) ) {
          continue; // implicit node reached
       }

       MEM = 1;
       // Taking Weiner links from current node to visit nodes at string depth++
       // TODO: Push the largest interval first to limit the size of the stack
       // TODO: use as alphabet the distinct symbols appearing in the range

       // last alphabet symbol regarded as separator
       for (ulint j=0;j<alphabet.size()-1; j++) { 
          node.first = tidx.left_extension(alphabet[j],sample); 
          node.second = qidx.left_extension(alphabet[j],qsample);   
          S.push(node);       
          dS.push(d+1);
          if (f) {
             fS.push(fidx.left_extension(alphabet[j],fsample));
          }
          
          // range not splitting, no MEM reported
          if (sample.size()+qsample.size()==node.first.size()+node.second.size()) 
             MEM = 0;
       }
       if (d > maxMEM)
          maxMEM = d;
       if (MEM and d>=kappa and output.is_open())
          if (!asymmetric)
             reportMEMs<T,TS>(tidx,qidx,sample,qsample,d,output,Bsuffix,Bprefix);
          else if (!f or fsample.is_invalid()) // MEM string not found in filter index, so not a dublicate
             reportAMEMs<T,TS>(tidx,qidx,sample,qsample,d,output);   
    }
    return maxMEM;
}

template<class T,class TS>
void find_mems(ifstream& in, ifstream& qin, ifstream& efg_in)
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    using std::chrono::microseconds;


    auto tstart = high_resolution_clock::now();
    string header;
    getline(efg_in,header); // header >nodes
    string nodes;
    getline(efg_in,nodes);
    getline(efg_in,header); // header >edges
    string edges;
    getline(efg_in,edges);
    getline(efg_in,header); // header >paths    
    string paths;
    getline(efg_in,paths);
    
    string nodes_without_gt(nodes);    
    nodes_without_gt.erase(std::remove(nodes_without_gt.begin(), nodes_without_gt.end(), '>'), nodes_without_gt.end());
    string edges_without_gt(edges);    
    edges_without_gt.erase(std::remove(edges_without_gt.begin(), edges_without_gt.end(), '>'), edges_without_gt.end());
    string paths_without_gt(paths);    
    paths_without_gt.erase(std::remove(paths_without_gt.begin(), paths_without_gt.end(), '>'), paths_without_gt.end());

    //bool* Be_suffix = new bool[edges_without_gt.size()];
    bit_vector Be_suffix(edges_without_gt.size(),0);
    //bool* Be_suffix_bwt = new bool[edges_without_gt.size()+1];    
    //bool* Be_prefix = new bool[edges_without_gt.size()];    
    bit_vector Be_prefix(edges_without_gt.size(),0);    
    //bool* Be_prefix_bwt = new bool[edges_without_gt.size()+1];        
    //bool* Bp_suffix = new bool[paths_without_gt.size()];
    bit_vector Bp_suffix(paths_without_gt.size(),0);
    //bool* Bp_suffix_bwt = new bool[paths_without_gt.size()+1];  
    //bool* Bp_prefix = new bool[paths_without_gt.size()];    
    bit_vector Bp_prefix(paths_without_gt.size(),0);    
    //bool* Bp_prefix_bwt = new bool[paths_without_gt.size()+1];    
    
    if (alphabet.size()==0) {
       alphabet = "ACGT#";
    }
    
    bool suffix = true;
    bool prefix = false;
    ulint j=0;
    for (ulint i=0; i<edges.size(); i++)
       if (edges[i]==alphabet[alphabet.size()-1]) { // edge boundary
          suffix = true; 
          prefix = false;
          Be_suffix[j]=suffix;
          Be_prefix[j++]=prefix;          
       } 
       else if (edges[i]!='>') {
          Be_suffix[j]=suffix;
          Be_prefix[j++]=prefix;          
       }
       else { // node boundary
          suffix = false;
          prefix = true;
       }
           
    suffix = true;
    prefix = false;
    j = 0;
    for (ulint i=0; i<paths.size(); i++)
       if (paths[i]==alphabet[alphabet.size()-1]) { // path boundary
          suffix = true; 
          prefix = false;
          Bp_suffix[j]=suffix;
          Bp_prefix[j++]=prefix;          
       } 
       else if (paths[i]!='>') {
          Bp_suffix[j]=suffix;
          Bp_prefix[j++]=prefix;          
       }
       else if (suffix)  // first node boundary
          suffix = false;
       else // second node boundary
          prefix = true;

    /* Converting bitvectors to BWT order, not used in the naive cross-product
    j = eidx.get_terminator_position();    
    j = eidx.LF(j);
    for (ulint i=0; i<edges_without_gt.size(); i++) {
       Be_suffix_bwt[j] = Be_suffix[nodes_without_gt.size()-i-1];
       j = eidx.LF(j);       
    } 
    delete[] Be_suffix;
    
    j = eidx.get_terminator_position(true);    
    j = eidx.LFR(j);
    for (ulint i=0; i<edges_without_gt.size(); i++) {
       Be_prefix_bwt[j] = Be_prefix[edges_without_gt.size()-i-1];
       j = eidx.LFR(j);       
    } 
    delete[] Be_prefix; 
   
    j = pidx.get_terminator_position();    
    j = pidx.LF(j);
    for (ulint i=0; i<paths_without_gt.size(); i++) {
       Bp_suffix_bwt[j] = Bp_suffix[paths_without_gt.size()-i-1];
       j = pidx.LF(j);       
    } 
    delete[] Bp_suffix;
    
    j = pidx.get_terminator_position(true);    
    j = pidx.LFR(j);
    for (ulint i=0; i<paths_without_gt.size(); i++) {
       Bp_prefix_bwt[j] = Bp_prefix[paths_without_gt.size()-i-1];
       j = pidx.LFR(j);       
    } 
    delete[] Bp_prefix;
    */
          
    // releasing memory
    nodes.clear();
    edges.clear();
    paths.clear();

    
    // Building indexes
    T nidx = T(nodes_without_gt);
    T eidx = T(edges_without_gt);
    T pidx = T(paths_without_gt);

    // releasing memory
    nodes_without_gt.clear();
    edges_without_gt.clear();
    paths_without_gt.clear();

    T idx;
    if (filter)
       idx.load(in);
    
    T qidx;
    qidx.load(qin);      

    ofstream output;
    if (output_file.size()!=0) {     
       output.open(output_file);
       if (!output) {
          cout << "Could not open output file " << output_file << endl;
       }
    }

    auto tmem = high_resolution_clock::now();       
    cout << "Exploring MEMs " << endl;

    ulint maxMEM;
    output << ">nodes" << endl;
    maxMEM = explore_mems<T,TS>(nidx,qidx,idx,output);
    cout << "Maximum node MEM is of length " << maxMEM << endl;
    output << ">edges" << endl;
    if (!asymmetric)
       maxMEM = explore_mems<T,TS>(eidx,qidx,idx,output,false,Be_suffix,Be_prefix);    
    else  // using node index as filter 
       maxMEM = explore_mems<T,TS>(eidx,qidx,nidx,output,true,Be_suffix,Be_prefix);    
    cout << "Maximum edge MEM is of length " << maxMEM << endl;
    //delete[] Be_suffix;    
    //delete[] Be_prefix;
    output << ">paths" << endl;
    if (!asymmetric) // using text index as filter, if filter=true
       maxMEM = explore_mems<T,TS>(pidx,qidx,idx,output,filter,Bp_suffix,Bp_prefix);        
    else // using edge index as filter
       maxMEM = explore_mems<T,TS>(pidx,qidx,eidx,output,true,Bp_suffix,Bp_prefix);           
    cout << "Maximum path MEM (spanning 3 nodes) is of length " << maxMEM << endl;
    //delete[] Bp_suffix;    
    //delete[] Bp_prefix;    
    // TODO
    // maxMEM = exploreLongPathMEMs();
    // cout << "Maximum path MEM (spanning more than 3 nodes) is of length " << maxMEM << endl;
       
    output.close();

    auto tend = high_resolution_clock::now();
    ulint tot_time = duration_cast<microseconds>(tend-tstart).count(); 
    ulint mem_time = duration_cast<microseconds>(tend-tmem).count(); 
    cout << "MEM exploration time     : " << mem_time << " microseconds" << endl;
    cout << "Total time     : " << tot_time << " microseconds" << endl;
}



int main(int argc, char** argv)
{
    if (argc < 3) help();

    int ptr = 1;

    while (ptr < argc - 3) 
       parse_args(argv, argc, ptr);

    string idx_file(argv[ptr]);
    string query_file(argv[ptr+1]);
    string efg(argv[ptr+2]);
    
    
    cout << "Loading and creating br-indexes" << endl;

    ifstream in;
    if (filter)
       in = ifstream(idx_file);
    ifstream qin(query_file);
    ifstream efg_in(efg);
       
    if (nplcp) {
       find_mems<br_index_nplcp<>,br_sample_nplcp >(in, qin, efg_in);
    }
    else {
       find_mems<br_index<>,br_sample >(in, qin, efg_in);
    }

    in.close();
    qin.close();
    efg_in.close();

}
