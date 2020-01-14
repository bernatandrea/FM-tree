
FMtree: A fast locating algorithm of FM-indexes for genomic data
============

## Authors
- [Robert Jambrečić](https://github.com/robijam)
- [Anel Hadžimuratagić](https://github.com/anelka19)
- [Andrea Bernat](https://github.com/bernatandrea)




Introduction
-------  

Here are the implementation of "FMtree: A fast locating algorithm of FM-indexes for genomic data".
This program was created as a project assignment for <a href="https://www.fer.unizg.hr/predmet/bio" target="_blank">Bioinformatics</a> class at Faculty of Electrical Engineering and Computing, University of Zagreb.







How to use
-------



For FM-tree:

* First to generate binary write this command in your command shell positioned inside the project.

Create executable file:
```
make 
```
Now you can run program:
```
FMtree.exe 
```

* After you run the program you will be asked to input text file you would like to search and sampling distance of suffix array SA.
`Please note that input text file cannot include any character which does not belong to {A, C, G, T, $(end)}`

* For generating the pattern you would like to search over your input text file choose the first option.
`Like the input text file, pattern file can only include the characters which belong to {A, C, G, T}. Note that there is no $ in pattern`. 
If some patterns consist of such characters, the results of FMtree would be incorrect.

* To build the index for the input text(genome) you will choose second option.
 For example, consider a text "ecoli.txt", its index consists of  BWT,SA,C(occ) arrays which are created by choosing 
 second option.

* For searching with FMtree choose third option. Result with count and locations will be saved in `results.txt` and 
program preferences will be written on a console.

* Example: After running FMtree it will report the following information:

     ![example1](https://github.com/bernatandrea/FM-tree/raw/master/example1.png) 


Test data:

* Test data are saved in folder `test`, were you can find input text tests from 100 - 1 000 000 characters 

* Ecoli genome also in `test/ecoli.txt` with more than 2 000 000 characters




Note
-------
* We adopt the divsufsort library[1] to build the suffix array, and build BWT from suffix array. As such when building the index, the memory requirement of FMtree is about 5 times larger than that of the input text.

* We also adopt wavalet tree algorithm[2] to use rank function which was causing the most time cost when it was linearly implemented.

## License

[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)

- **[MIT license](http://opensource.org/licenses/mit-license.php)**


References
-------
[1] [https://github.com/felipelouza/bwt-lcp-in-place/tree/master/external/divsufsort-lite](https://github.com/felipelouza/bwt-lcp-in-place/tree/master/external/divsufsort-lite)
[2][https://n1chre.github.io/bioinf/](https://n1chre.github.io/bioinf/)
