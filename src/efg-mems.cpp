#include <iostream>
#include <chrono>
#include <cstdlib>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include <sdsl/bit_vectors.hpp>
#include "br_index.hpp"
#include "utils.hpp"


using namespace bri;
using namespace std;
using namespace sdsl;

// Global variables assigned with optional parameters
ulint kappa = 1; // MEM threshold
bool asymmetric = false;
bool filter = false;
string alphabet = "";
string output_file = "";

void help()
{
	cout << "efg-mems: locate all MEMs between queries and elastic founder graph" << endl;

	cout << "Usage: efg-mems [options] <text> <queries> <efg>" << endl;
    cout << "   --asymmetric  use asymmetric definition for MEMs " << endl;
	cout << "   --filter      filter out MEMs that do not occur in <text> " << endl;
	cout << "   -k      MEM threshold" << endl;
	cout << "   -a      alphabet string with last three symbols regarded special, default ACGTN#0" << endl;	
	cout << "   -o      output file where lines x,i,d are outputed for MEMs Q[i..i+d-1]=T[x..x+d-1]; " << endl;
    cout << "           in symmetric mode (default) the matches are locally maximal, i.e., Q[i-1]!=T[x-1] and Q[i+d]!=T[x+d] " << endl;
	cout << "           in asymmetric mode the matches are globally maximal, i.e., Q[i..i+d] or Q[i-1..i+d-1] do not " << endl;
	cout << "           occur in T; only one occurrence T[x..x+d-1] is reported in asymmetric mode" << endl;		
	cout << "           output is formatted as fasta file with MEMs listed according to the input efg fasta" << endl;		
	cout << "   -f     file containing alphabet " << endl;
	cout << "   <text>  index file of the concatenation of MSA rows (with extension .bri)" << endl;
	cout << "   <queries>  index file of concatenation of queries (with extension .bri)" << endl;
	cout << "               concatenation format: #AGGATG#AGATGT#, where # is separator symbol" << endl;		
	cout << "   <efg>      elastic founder graph in GFA format (with extension .gfa)" << endl;
	exit(0);
}

bool parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;


    if (s.compare("--asymmetric") == 0)
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

void read_gfa(ifstream& gfa, string& nodes, string& edges, string& paths) {

    vector<string> node_labels;
    vector<string> edge_labels;
    vector<int> edge_node_start;
    vector<int> edge_node_end;    
    vector<string> path_labels;
    unordered_map<int, unordered_set<char>> lext; // chars to the left from nodes 
    unordered_map<int, unordered_set<char>> rext;  // chars to the right from nodes
    unordered_map<int, unordered_set<int>> adj_list;
    
    string lextnodes, rextnodes, lextedges, rextedges;
    
    string line;
    while (getline(gfa, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip empty or comment lines
        }

        vector<string> fields;
        string field;
        for (char c : line) {
            if (c == '\t') {
                fields.push_back(field);
                field.clear();
            } else {
                field += c;
            }
        }
        fields.push_back(field);

        if (fields[0] == "S") {
            node_labels.push_back(fields[2]);
        } else if (fields[0] == "L") {
            int start_node_id = stoi(fields[1]);
            bool start_node_orientation = (fields[2] == "-");
            int end_node_id = stoi(fields[3]);
            bool end_node_orientation = (fields[4] == "-");
            lext[end_node_id].insert(node_labels[start_node_id-1][node_labels[start_node_id-1].size()-1]);
            rext[start_node_id].insert(node_labels[end_node_id-1][0]);            
            adj_list[start_node_id].insert(end_node_id);
            if (!end_node_orientation) {
                adj_list[end_node_id].insert(start_node_id);
            }
            string start_node_label = node_labels[start_node_id - 1];
            if (start_node_orientation) {
                start_node_label = start_node_label + "_rev";
            }
            string end_node_label = node_labels[end_node_id - 1];
            if (end_node_orientation) {
                end_node_label = end_node_label + "_rev";
            }
            edge_labels.push_back(start_node_label + ">" +end_node_label);
            edge_node_start.push_back(start_node_id);
            edge_node_end.push_back(end_node_id);
        }
    }

    nodes = alphabet[alphabet.size()-1];
    lextnodes = alphabet[alphabet.size()-1];
    rextnodes = alphabet[alphabet.size()-1];
    for (int i=0; i < node_labels.size(); i++) {
       nodes += alphabet[alphabet.size()-1]+node_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // left and right chars not yet known
       if (i==0) 
          lextnodes += alphabet[alphabet.size()-3] + node_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // ambiguous left
       else if (lext[i].size()==1) { 
          for (char lextchar : lext[i])
             lextnodes += lextchar + node_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // unique left
       }
       else
          lextnodes += alphabet[alphabet.size()-3] + node_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // ambiguous left
       if (rext[i].size()==1) {
          for (char rextchar : rext[i])
             rextnodes += alphabet[alphabet.size()-1]  + node_labels[i] + rextchar + alphabet[alphabet.size()-1]; // unique right
       }
       else
          rextnodes += alphabet[alphabet.size()-1] + node_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; //ambiguous right
    }
  
    edges = alphabet[alphabet.size()-1];
    lextedges = alphabet[alphabet.size()-1];
    rextedges = alphabet[alphabet.size()-1];
    for (int i=0; i < edge_labels.size(); i++) {
       edges += alphabet[alphabet.size()-1]+edge_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // left and right chars not yet known
       if (i==0) 
          lextedges += alphabet[alphabet.size()-3] + edge_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // ambiguous left
       else if (lext[edge_node_start[i]].size()==1) { 
          for (char lextchar : lext[edge_node_start[i]])
             lextedges += lextchar + edge_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // unique left
       }
       else
          lextedges += alphabet[alphabet.size()-3] + edge_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; // ambiguous left
       if (rext[edge_node_end[i]].size()==1) {
          for (char rextchar : rext[edge_node_end[i]])
             rextedges += alphabet[alphabet.size()-1]  + edge_labels[i] + rextchar + alphabet[alphabet.size()-1]; // unique right
       }
       else
          rextedges += alphabet[alphabet.size()-1] + edge_labels[i] + alphabet[alphabet.size()-1] + alphabet[alphabet.size()-1]; //ambiguous right
    }
    int minlen = 1000;   
    for (auto& kv1 : adj_list) {
        int node1_id = kv1.first;
        for (int node2_id : kv1.second) {
            for (int node3_id : adj_list[node2_id]) {
                string node1_label = node_labels[node1_id - 1];
                string node2_label = node_labels[node2_id - 1];
                string node3_label = node_labels[node3_id - 1];
                if (node1_label.size()+node2_label.size()+node2_label.size() < minlen)
                   minlen = node1_label.size()+node2_label.size()+node2_label.size();
                path_labels.push_back(node1_label + ">" + node2_label + 
                                                    ">" + node3_label);
            }
        }
   }

   cout << "minlen " << minlen << endl;
   paths = alphabet[alphabet.size()-1];
   for (string path_label : path_labels) {
       paths += alphabet[alphabet.size()-1] + path_label + alphabet[alphabet.size()-1] +alphabet[alphabet.size()-1]; // ambiguous left and right
   }
   // combining left and right chars to nodes and edges
   nodes[1] = lextnodes[1];
   nodes[nodes.size()-2] = rextnodes[nodes.size()-2];
   for (int i=2; i<nodes.size()-2; i++)
      if (nodes[i+1]==alphabet[alphabet.size()-1]) {
         nodes[i] = rextnodes[i];
         nodes[i+2] = lextnodes[i+2];
      }
   edges[1] = lextedges[1];
   edges[edges.size()-2] = rextedges[edges.size()-2];
   for (int i=2; i<edges.size()-2; i++)
      if (edges[i+1]==alphabet[alphabet.size()-1]) {
         edges[i] = rextedges[i];
         edges[i+2] = lextedges[i+2];
      }
}


template <class T, class TS>
bool is_right_maximal(T idx, TS sample) 
{
   uchar c = idx.bwt_at(sample.rangeR.first,true);
   TS right = idx.right_extension(c,sample);
   if (sample.size() == right.size())
      return 0;
   else
      return 1;
}
    
template <class T, class TS>    
bool is_left_maximal(T idx, TS sample) 
{
   uchar c = idx.bwt_at(sample.rangeR.first,false);
   TS left = idx.left_extension(c,sample);
   if (sample.size() == left.size())
      return 0;
   else
      return 1;
}

template<class T, class TS>
void reportMEMs(T idx, T qidx, TS sample, TS qsample, ulint d, ofstream& output, bit_vector Bsuffix = bit_vector(), bit_vector Bprefix = bit_vector(), ulint* d_bwt = NULL, ulint* saidx = NULL)
{ 
   // a bit naive implementation of the cross product 
   // additive factor alphabet.size()^2(t_biBWTstep+\alphabet.size()^2) slower than an optimal implementation
   std::vector<ulint>** a= new  std::vector<ulint>*[alphabet.size()];
   std::vector<ulint>** b= new  std::vector<ulint>*[alphabet.size()];
   TS** Sa= new  TS*[alphabet.size()];
   TS** Sb= new  TS*[alphabet.size()];
   bool** Ia= new  bool*[alphabet.size()];
   bool** Ib= new  bool*[alphabet.size()];

   TS left,right;
   for (ulint i=0; i<alphabet.size()-1; i++) { // not branching on graph concat separator char, as no valid match can start/end at the boundaries
      a[i] = new std::vector<ulint>[alphabet.size()];
      b[i] = new std::vector<ulint>[alphabet.size()];      
      Sa[i] = new TS[alphabet.size()];
      Sb[i] = new TS[alphabet.size()];
      Ia[i] = new bool[alphabet.size()];
      Ib[i] = new bool[alphabet.size()];      
      for (ulint j=0; j<alphabet.size()-1; j++) {
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
   for (ulint i=0; i<alphabet.size()-1; i++) 
      for (ulint j=0; j<alphabet.size()-1; j++)
         if (Ia[i][j]) {
            for (ulint ii=0; ii<alphabet.size()-1; ii++)
               for (ulint jj=0; jj<alphabet.size()-1; jj++)
                  if (i!=ii and j!=jj and Ib[ii][jj]) {
                     if (a[i][j].size()==0)  // locate charged on the output
                        a[i][j] = qidx.locate_sample(Sa[i][j]);
                     if (b[ii][jj].size()==0)  // locate charged on the output
                        if (saidx==NULL) // no suffix array, using slower locate
                           b[ii][jj] = idx.locate_sample(Sb[ii][jj]);
                        else 
                           for (ulint jjj=0; jjj<Sb[ii][jj].range.second-Sb[ii][jj].range.first+1; jjj++)
                              b[ii][jj].push_back(saidx[Sb[ii][jj].range.first+jjj]);              
                     for (ulint iii=0; iii<a[i][j].size(); iii++)  {
                        for (ulint jjj=0; jjj<b[ii][jj].size(); jjj++)
                           //if (d_bwt!=NULL and d_bwt[Sb[ii][jj].range.first+jjj]<=d+1)
                           //   output << saidx[Sb[ii][jj].range.first+jjj]+1 << "," << a[i][j][iii]+1 << "," << d << endl;
                           //else 
                           //   output << b[ii][jj][jjj]+1 << "," << a[i][j][iii]+1 << "," << d << endl;                           
                              
                           if (Bsuffix.size()==0 or Bprefix.size()==0 or (Bsuffix[b[ii][jj][jjj]+1] and Bprefix[b[ii][jj][jjj]+d])) {
                              output << b[ii][jj][jjj]+1 << "," << a[i][j][iii]+1 << "," << d << endl;
                           }
                     }
                  }   
         }

   for (ulint i=0; i<alphabet.size()-1; i++) {
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
   for (ulint i=0; i<alphabet.size()-1; i++) {
      a[i] = new std::vector<ulint>[alphabet.size()];
      b[i] = new bool[alphabet.size()];
      Sa[i] = new TS[alphabet.size()];
      Ia[i] = new bool[alphabet.size()];
      for (ulint j=0; j<alphabet.size()-1; j++) {
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
   for (ulint i=0; i<alphabet.size()-1; i++) 
      for (ulint j=0; j<alphabet.size()-1; j++)
         if (Ia[i][j] and b[i][j]) { 
            a[i][j] = qidx.locate_sample(Sa[i][j]);
            for (ulint iii=0; iii<a[i][j].size(); iii++) 
               output << a[i][j][iii]+1 << ","  << sample.j -sample.d << "," << d << endl;          
         }
   for (ulint i=0; i<alphabet.size()-1; i++) {
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
ulint explore_mems(T tidx, T qidx, T fidx, ofstream& output, bool f = false, bit_vector Bsuffix = bit_vector(), bit_vector Bprefix = bit_vector(), ulint* d_bwt = NULL, ulint* sa = NULL)
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
       if ((!is_right_maximal<T,TS>(tidx,sample) and !is_right_maximal<T,TS>(qidx,qsample)) and 
           tidx.bwt_at(sample.rangeR.first,true)==qidx.bwt_at(qsample.rangeR.first,true) ) {
          continue; // implicit node reached
       }

       MEM = 1;
       // Taking Weiner links from current node to visit nodes at string depth++
       // TODO: Push the largest interval first to limit the size of the stack
       // TODO: use as alphabet the distinct symbols appearing in the range

       // not branching with three last alphabet characters
       for (ulint j=0;j<alphabet.size()-3; j++) { 
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
             reportMEMs<T,TS>(tidx,qidx,sample,qsample,d,output,Bsuffix,Bprefix,d_bwt,sa);
          else if (!f or fsample.is_invalid()) // MEM string not found in filter index, so not a dublicate
             reportAMEMs<T,TS>(tidx,qidx,sample,qsample,d,output);   
    }
    return maxMEM;
}

template<class T,class TS>
void find_mems(ifstream& in, ifstream& qin, ifstream& efg_in,string path_prefix)
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    using std::chrono::microseconds;


    auto tstart = high_resolution_clock::now();
    
    string nodes, edges, paths;

    if (alphabet.size()==0) {
       alphabet = "ACGTN#$";
       // last char is a separator in graph concatenation, omitted in all algorithms
       // second last char is separator in queries, not used in MEM exploration, but used in cross-product
       // third last char marks ambiguous left- and right-extension, not used in MEM exploration, but used in cross-product
    }

    read_gfa(efg_in,nodes,edges,paths);
         
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
           
    ulint* d_edge = new ulint[edges_without_gt.size()+1]; // distance from pos in start node to the start of end node
    j = 0;
    for (ulint i=edges_without_gt.size()-1; i!=0; i--) {
       if (Be_prefix[i])
          j = 0;
       else 
          j++;
       if (!Be_suffix[i])
          d_edge[i]=edges_without_gt.size();
       else
          d_edge[i]=j;
    }
    d_edge[edges_without_gt.size()]=edges_without_gt.size();
    
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

    ulint* d_path = new ulint[paths_without_gt.size()+1]; // distance from pos in start node to the start of end node
    j = 0;
    for (ulint i=paths_without_gt.size()-1; i!=0; i--) {
       if (Bp_prefix[i])
          j = 0;
       else 
          j++;
       if (!Bp_suffix[i])    
          d_path[i]=paths_without_gt.size();
       else
          d_path[i]=j;       
    }
    d_path[paths_without_gt.size()]=paths_without_gt.size();

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

    
    // Building and saving indexes if they don't exist
    ifstream nidx_in(path_prefix+".nodes.bri");
    T nidx;
    if (nidx_in.is_open())
       nidx.load(nidx_in);
    else {
       nidx = T(nodes_without_gt);
       nidx.save_to_file(path_prefix+ ".nodes"); 
    }
    ifstream eidx_in(path_prefix+".edges.bri");
    T eidx;
    if (eidx_in.is_open())
       eidx.load(eidx_in);
    else {
       eidx = T(edges_without_gt);
       eidx.save_to_file(path_prefix+".edges");
    } 
    ifstream pidx_in(path_prefix+".paths.bri");
    T pidx;
    if (pidx_in.is_open())
       pidx.load(pidx_in);
    else {
       pidx = T(paths_without_gt);
       pidx.save_to_file(path_prefix+".paths"); 
    }
    
    /* Converting d_edge and d_path to BWT order */
    // Computing suffix arrays also
    ulint* d_edge_bwt = new ulint[edges_without_gt.size()+1];    
    ulint* sa_edge = new ulint[edges_without_gt.size()+1];
    j = eidx.get_terminator_position();   
    d_edge_bwt[j]=d_edge[0];
    sa_edge[j]=0;
    j = eidx.LF(j);
    for (ulint i=1; i<=edges_without_gt.size(); i++) {
       d_edge_bwt[j] = d_edge[edges_without_gt.size()-i+1];
       sa_edge[j]=edges_without_gt.size()-i+1;
       j = eidx.LF(j);    
    } 
    delete[] d_edge;
    
    ulint* d_path_bwt = new ulint[paths_without_gt.size()+1];    
    ulint* sa_path = new ulint[paths_without_gt.size()+1];    
    j = pidx.get_terminator_position(); 
    d_path_bwt[j]=d_path[0]; 
    sa_path[j]=0;
    j = pidx.LF(j);
    for (ulint i=1; i<=paths_without_gt.size(); i++) {
       d_path_bwt[j] = d_path[paths_without_gt.size()-i+1];
       sa_path[j]=paths_without_gt.size()-i+1;
       j = pidx.LF(j);       
    } 
    delete[] d_path;
    
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
       maxMEM = explore_mems<T,TS>(eidx,qidx,idx,output,false,Be_suffix,Be_prefix,d_edge_bwt,sa_edge);    
    else  // using node index as filter 
       maxMEM = explore_mems<T,TS>(eidx,qidx,nidx,output,true,Be_suffix,Be_prefix,d_edge_bwt,sa_edge);    
    cout << "Maximum edge MEM is of length " << maxMEM << endl;
    //delete[] Be_suffix;    
    //delete[] Be_prefix;
    output << ">paths" << endl;
    string lextpaths = "";
    string rextpaths = "";
    if (!asymmetric) // using text index as filter, if filter=true
       maxMEM = explore_mems<T,TS>(pidx,qidx,idx,output,filter,Bp_suffix,Bp_prefix,d_path_bwt,sa_path);        
    else // using edge index as filter
       maxMEM = explore_mems<T,TS>(pidx,qidx,eidx,output,true,Bp_suffix,Bp_prefix,d_path_bwt,sa_path);           
    cout << "Maximum path MEM (spanning 3 nodes) is of length " << maxMEM << endl;
    //delete[] Bp_suffix;    
    //delete[] Bp_prefix;    
    // TODO
    // maxMEM = exploreLongPathMEMs();
    // cout << "Maximum path MEM (spanning more than 3 nodes) is of length " << maxMEM << endl;
       
    output.close();
    
    delete[] d_edge_bwt;
    delete[] d_path_bwt;
    delete[] sa_edge;
    delete[] sa_path;
    
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
       
    find_mems<br_index<>,br_sample >(in, qin, efg_in, efg);

    in.close();
    qin.close();
    efg_in.close();

}
