https://travis-ci.org/PollockLab/pclouds.svg?branch=feature%2Ftravis-ci

PClouds
-------

Goal: Identify repeat structure in large eukaryotic genomes using oligonucleotide counts. It works efficiently on a single desktop computer with 1 Gb memory.

Abstract
--------

The P-clouds package is designed to identify repeat structure in large eukaryotic genomes using oligonucleotide counts. It works efficiently on a single desktop computer with 1 Gb memory. The basic program is described in Gu et al. (2008), below, with more details, analysis of human genome, and description of element-specific P-clouds in de Koning et al. (2011).

Compatibility
-------------

The program should should be compatible with any system that has the following:

-   git
-   cmake
-   C++11 compatible C++ compiler

How to use
----------

First you use cmake to create the make files:

    cd build
    cmake ..

then you compile it:

    make

and then change the control file to your needs, and run it:

    ./pclouds

