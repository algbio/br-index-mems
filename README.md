# br-index

## News

New features in this fork:

- Support for MEM finding 

## About

This repository provides the __bi-directional__ r-index (_br-index_).

The r-index is the compressed text index which supports efficient count(P) and locate(P).
Its uses O(r) words of space, whose r is the number of equal-letter runs in BWT of the text.

The br-index we have proposed is achieved by using the mechanism of the r-index but adding new structures. It allows for bi-directional extension of patterns, _left-extension_ and _right-extension_. Also, you can locate all the occurrences of the current pattern at any step of the search.

## System Requirements

- This project is based on [sdsl-lite](https://github.com/simongog/sdsl-lite) library.
Install sdsl-lite beforehand and modify variables SDSL_INCLUDE and SDSL_LIB in _CMakeLists.txt_.

- This project has been tested under gcc 4.8.5 and gcc 7.5.0.

## How to Use

Firstly, clone the repository. Since a submodule is used ([iutest](https://github.com/srz-zumix/iutest)), recursive cloning is necessary.
```bash
git clone --recursive https://github.com/U-Ar/br-index.git
```
In order to build, execute following commands: (This project is using CMake)
```bash
mkdir build
cd build
cmake ..
make
```
7 executables will be created in the _build_ directory.
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
	<dt>run_tests</dt>
	<dd>runs unit tests.</dd>
</dl>

Also, you can run unit tests by
```bash
make test-bri
```

## Versions

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

Cite the following paper:
- Arakawa, Y., Navarro, G., & Sadakane, K. (2022). Bi-Directional r-Indexes. In 33rd Annual Symposium on Combinatorial Pattern Matching (CPM 2022). Schloss Dagstuhl-Leibniz-Zentrum für Informatik.

It is more desirable to cite the following papers in addition, which are the original papers of the r-index:
- Gagie, T., Navarro, G., & Prezza, N. (2018). Optimal-time text indexing in BWT-runs bounded space. In Proceedings of the Twenty-Ninth Annual ACM-SIAM Symposium on Discrete Algorithms (pp. 1459-1477). Society for Industrial and Applied Mathematics.
- Gagie, T., Navarro, G., & Prezza, N. (2020). Fully functional suffix trees and optimal text searching in BWT-runs bounded space. Journal of the ACM (JACM), 67(1), 1-54.
