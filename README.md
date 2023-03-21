# br-index-mems

## News

This is fork of br-index by Yuma Arakawa, adding support for MEM finding:

- Support for finding maximal exact matches (MEMs) between a set of query patterns and a text
- Support for finding MEMS between a set of query patterns and an elastic founder graph


## About the root repository

The root repository provides the __bi-directional__ r-index (_br-index_).

The r-index is the compressed text index which supports efficient count(P) and locate(P).
Its uses O(r) words of space, where r is the number of equal-letter runs in BWT of the text.

The br-index is achieved by using the mechanism of the r-index but adding new structures. It allows for bi-directional extension of patterns, _left-extension_ and _right-extension_. Also, you can locate all the occurrences of the current pattern at any step of the search.

## System Requirements

- This project is based on [sdsl-lite](https://github.com/simongog/sdsl-lite) library.
Install sdsl-lite beforehand and modify variables SDSL_INCLUDE and SDSL_LIB in _CMakeLists.txt_.

- This project has been tested under gcc 4.8.5 and gcc 7.5.0.

## How to Use

Firstly, clone the repository. Since a submodule is used ([iutest](https://github.com/srz-zumix/iutest)), recursive cloning is necessary.
```bash
git clone --recursive https://github.com/velimakinen/br-index-mems.git
```
In order to build, execute following commands: (This project is using CMake)
```bash
mkdir build
cd build
cmake ..
make
```
8 executables will be created in the _build_ directory.
<dl>
	<dt>bri-build</dt>
	<dd>Builds the br-index on the input text file.</dd>
	<dt>bri-locate</dt>
	<dd>Locates the occurrences of the given pattern using the index. Provide a pattern file in 
	the <a href="http://pizzachili.dcc.uchile.cl/experiments.html">Pizza&Chili format</a>. You can give an option "-m (number)" for the number of mismatched characters allowed (0 by default).</dd>
	<dt>bri-count</dt>
	<dd>Counts the number of the occurrences of the given pattern using the index. Its usage is same as bri-locate.</dd>
	<dt>bri-seedex</dt>
	<dd>Applies the seed-and-extend approach to the given pattern. Exactly matches the core region and extends with some mismatches.</dd>
	<dt>bri-space</dt>
	<dd>Shows the statistics of the text and the breakdown of the index space usage.</dd>
	<dt>bri-mem</dt>
	<dd>Finds MEMs between a collection of patterns and a text.</dd>
	<dt>bri-mem-efg</dt>
	<dd>Finds MEMs between a collection of patterns and an elastic founder graph.</dd>
	<dt>run_tests</dt>
	<dd>runs unit tests.</dd>
</dl>

You can run unit tests by
```bash
make test-bri
```

For MEM finding, you can test the following:
```bash
mkdir inputs
mkdir outputs
cp ../patterns.txt inputs/
cp ../text.txt inputs/
./br-build inputs/patterns.txt
./br-build inputs/text.txt
./bri-mem -k 4 -o outputs/MEMs.txt inputs/text.txt.bri inputs/patterns.txt.bri
cat outputs/MEMs.txt
```
<dt>This should print MEMs x,i,d s.t. P[i..i+d-1]=T[x..x+d-1] and d>=4.</dt>
<dt>In our case (note 0-based indexing):</dt>
<dt>15,2,6<dt>
<dt>65,2,6</dt>
<dt>31,18,4</dt>
<dt>25,1,10</dt>
<dt>12,3,4</dt>
<dt>60,1,6</dt>

## Versions (from root repository)

<dl>
	<dt>br_index.hpp (default)</dt>
	<dd>The simple implementation of br-index used in the experiments. Only the variables necessary for locate <i>(j,d,len)</i> are maintained, which are sufficient to compute <i>locate.</i></dd>
	<dt>br_index_nplcp.hpp</dt>
	<dd>The implementation without PLCP. <i>(p,j,d,len)</i> are maintained. It computes <i>locate</i> by calculating <i>p'=LF^d(p)</i> and comparing <i>p'</i> with
	[<i>s, e</i>]. Larger than the normal one while computing <i>locate</i> is faster when occ is large compared to |P|.
	<dt>br_index_naive.hpp (old)</dt>
	<dd>The naive implementation of br-index. All the variables <i>p,j,d,pR,jR,dR,len</i> are maintained during the search. Not space-efficient, implemented mainly for the educational purpose and the possible future use. (It's not updated now, so it doesn't function)</dd>
</dl>

## Citation

Cite the following paper for br-index:
- Arakawa, Y., Navarro, G., & Sadakane, K. (2022). Bi-Directional r-Indexes. In 33rd Annual Symposium on Combinatorial Pattern Matching (CPM 2022). Schloss Dagstuhl-Leibniz-Zentrum für Informatik.

It is more desirable to cite the following papers in addition, which are the original papers of the r-index:
- Gagie, T., Navarro, G., & Prezza, N. (2018). Optimal-time text indexing in BWT-runs bounded space. In Proceedings of the Twenty-Ninth Annual ACM-SIAM Symposium on Discrete Algorithms (pp. 1459-1477). Society for Industrial and Applied Mathematics.
- Gagie, T., Navarro, G., & Prezza, N. (2020). Fully functional suffix trees and optimal text searching in BWT-runs bounded space. Journal of the ACM (JACM), 67(1), 1-54.

MEM finding is analogous to Algorithm 11.3 at page 226 in 
- V. Mäkinen, D. Belazzougui, F. Cunial, A. I. Tomescu: Genome-Scale Algorithm Design. Cambridge University Press, 2015.
- See also: D. Belazzougui, F. Cunial, J. Kärkkäinen, V. Mäkinen:
Linear-time String Indexing and Analysis in Small Space. ACM Trans. Algorithms 16(2): 17:1-17:54 (2020)

MEM finding on EFGs is analogous to Theorem 13 in 
- Nicola Rizzo, Manuel Cáceres, Veli Mäkinen: Chaining of Maximal Exact Matches in Graphs. https://arxiv.org/abs/2302.01748

