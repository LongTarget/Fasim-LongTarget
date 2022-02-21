# Fasim-LongTarget: A tool for genome-wide accurate prediction of lncRNA/DNA binding

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
    + [Compilation](#compilation)
    + [Running](#running)
    + [Help information](#help-information)
    + [Time consumption](#time-consumption)
- [Demo](#demo)
    + [Inputs and their formats](#inputs-and-their-formats)
    + [Results](#results)
    + [Example datasets](#example-datasets)
- [Bug reports](#bug-reports)
- [License](./LICENSE)
- [Citation](#citation)

# Overview
Many lncRNAs can bind to DNA sequences by forming RNA:DNA triplexes, by which they recruit histone modification and DNA methylation enzymes to binding sites to modify the epigenome. Since the formation of triplexes between lncRNA and DNA sequences follows Hoogsteen and reverse Hoogsteeen base-pairing rules, lncRNA/DNA bindings can be predicted computationally.

LongTarget has been developed to predict one or many lncRNA's DNA binding motifs and binding sites in one or many genome regions based on all known ncRNA/DNA base pairing rules. The fasim-LongTarget is a variant of LongTarget that was modified to reduce the time consumption. Fasim-LongTarget consists of a few C/C++ programs, and is distributed under the AGPLv3 license. Tests on experimentally-generated lncRNA/DNA binding datasets indicate the good performance and low time consumption of fasim-LongTarget.

# Repo Contents
- [fasim-LongTarget.cpp](./fasim-LongTarget.cpp): The main program for generating the executable "fasim".
- [rules.h](./rules.h): Base-pairing rules and codes for handling these rules.
- [sim.h](./sim.h): The SIM program for local alignment.
- [stats.h](./stats.h): Michael Farrar's code (with SSE2) for local alignment.
- [sswNew.cpp](./sswNew.cpp): The Program for getting the maximal score of each column, where the maximal scores are > threshold.
- [ssw.h](./ssw.h): The header file of sswNew.cpp.
- [ssw_cpp.cpp](./ssw_cpp.cpp): The program for getting the maximal scores of local maximum.
- [ssw_cpp.h](./ssw_cpp.h): The header file of ssw_cpp.cpp.
- [fastsim.h](./fastsim.h): The fasim program for getting k non-overlapping local alignments.
- [H19.fa](./H19.fa): An example of lncRNA sequence.  
- [testDNA.fa](./testDNA.fa): An example of DNA sequence. 

# System Requirements
- OS: Linux. Fasim-LongTarget has been tested to compile and run under CentOS 6.0. 
- System software: g++.
- RAM: 16G or above, depending on the number of lncRNAs and length of genome region.
- CPU: 4 cores or above, depending on the number of lncRNAs and length of genome region.

# Installation Guide
## Compilation
Typically, this command will generate an executable fasim-LongTarget program (named fasim): 

```
g++ fasim-LongTarget.cpp ssw_cpp.cpp sswNew.cpp -O -msse2 -o fasim
```

## Running 
A simple case is:

```
./fasim -f1 testDNA.fa -f2 H19.fa -O output/ -lg 40 
```

A more complex case is:

```
./fasim -f testDNA.fa -s H19.fa -O output/ -c 6000 -i 70 -S 1.0 -nt 25 -na 1000 -pc 1 -pt -500 -ds 10 -lg 60
```

In this command, output path is output/, cut sequence's length is 6000, identity is 70%, stability is 1.0, ntMin is 25 nt, ntMax is 1000 nt, penaltyC is 1, penaltyT is -500, distance between TFOs is 10, min length of triplexes is 60 (for more details about parameters, visit the website http://lncRNA.smu.edu.cn).

## Help information
Here is a brief explanation of the command line arguments:

```
Options   Parameters      Functions
f1   DNA sequence file  A string, indicate the DNA sequence file name.
f2   RNA sequence file  A string, indicate the RNA sequence file name.
O    Output path        A string, indicate the directory into which the results are outputted.
c    Cutlength          An integer, indicate the length of each segment, the default value is 5000.
i    identity           An integer, indicate the criterion of alignment output, the default value is 60.
S    stability          A floating point, indicate the criterion of base-pairing, the default value is 1.0.
ni   ntmin              An integer, indicate the min length of triplexes, the default value is 20.
na   ntmax              An integer, indicate the max length of triplexes, the default is 100000 but is rarely used.
pc   penaltyC           An integer, indicate penalty, the default value is 0.
pt   penaltyT           An integer, indicate penalty, the default value is -1000.
ds   c_dd               An integer, indicate the distance between TFOs, the default value is 15.
lg   c_length           An integer, indicate the min length of triplexes, the default value is 50.
```

## Time consumption
This depends on the number and length of lncRNAs and the length of genome regions. The expected running time for a 3000bp lncRNA and a 5000bp DNA is ~2 seconds on a normal desktop computer. 

# Demo
## Inputs and their formats
H19.fa and testDNA.fa are two demo examples of lncRNA and DNA sequence files, respectively. To obtain more details, go to our website http://lncRNA.smu.edu.cn and/or check files in the "examples" subdirectory.

In the lncRNA sequence file, the title line should be like ">lncRNA_name", and the lncRNA sequence should be in a new line.

In the DNA sequence file, the title line should be like ">species|chr|start-end", where chr, start and end indicate the genomic coordinate, and the DNA sequence should be in a new line.

## Results
The results include three files whose filenames ending with: (1)*TFOsorted, (2)*TFOclass1, (3)*TFOclass2. The TFOsorted file contains the details of all triplexes, the TFOclass1 file contains the TTS distribution of TFO1 in the genome region, and the TFOclass2 file contains the TTS distribution of TFO2. 

## Example datasets
The computationally-predicted and experimentally-generated lncRNA/DNA binding datasets of MEG3, NEAT1 and MALAT1 are given in the subdirectory "examples".

# Bug reports
Please send comments and bug reports to: zhuhao@smu.edu.cn.

# License
The program is distributed under the AGPLv3 license.

# Citation
