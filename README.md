# Welcome to Optimised CVTree

### Description
Composition Vector Tree, or CVTree, is a useful tool in determining the phylogenetic relationships between genome sequences i.e., determining how similar genome sequences are to each other. A poor performance implementation of the CVTree software was safely parallelised and optimised for significant performance improvements (about 4.4 times speed-up ratio).

On top of parallelisation, this project focusses on utilisation of programming techniques to optimise cache access. This was achieved via refactoring code as well as employing suitable data structures where neccessary( e.g. vector of pairs).

Please read the attached report for a more detailed explaination.

* `original.cpp` contains the original poor performance sequential implementation  
* `improved.cpp` contains the parallelised and optimised program


### Learnings
    * Profiling
    * Von Neumann Architecture and Bottleneck concepts
    * Cache access optimisation
    * Analysing data dependencies
    * Data structures for cache optimisation
    * High performance computing
    * Amdahl's Law

### Technologies
    * OpenMP
    * C++
    * VS Profiler

## Compiling
Use MicroSoft Visual C++ (MSVC) compiler with Maximum Speed optimiser i.e., /O2

## Running
Run the program on Visual Studio after creating a project for it. 
Inside `improved.cpp`, there are two relevant functions `CompareAllBacteriaPar` and `CompareAllBacteriaSeq` indicating the parallel and sequential function respectively (by default `CompareAllBacteriaPar` is running).

## Hardware Requirements
This program exploits multiple CPU cores. It was originally tested on and configured to run on a system with a minimum of 12 Virtual Cores (6 Physical Cores) and 16 GB RAM.
