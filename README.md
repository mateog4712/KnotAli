# KnotAli

#### Description:
Software implementation of KnotAli.      
KnotAli is an algorithm for predicting the pseudoknotted secondary structures of RNA using relaxed Hierarchical Folding.

#### Supported OS: 
Linux, macOS


### Installation:  
Requirements: A compiler that supports C++11 standard (tested with g++ version 4.9.0 or higher), Pthreads, and CMake version 3.1 or greater.    

[CMake](https://cmake.org/install/) version 3.1 or greater must be installed in a way that HFold can find it.    
To test if your Mac or Linux system already has CMake, you can type into a terminal:      
```
cmake --version
```
If it does not print a cmake version greater than or equal to 3.1, you will have to install CMake depending on your operating system.

[Linux instructions source](https://geeksww.com/tutorials/operating_systems/linux/installation/downloading_compiling_and_installing_cmake_on_linux.php)

#### Steps for installation   
1. [Download the repository](https://github.com/mateog4712/KnotAli.git) and extract the files onto your system.
2. From a command line in the root directory (where this README.md is) run
```
cmake -H. -Bbuild
cmake --build build
```   
If you need to specify a specific compiler, such as g++, you can instead run something like   
```
cmake -H. -Bbuild -DCMAKE_CXX_COMPILER=g++
cmake --build build
```   
This can be useful if you are getting errors about your compiler not having C++11 features.

After installing you can move the executables wherever you wish, but you should not delete or move the simfold folder, or you must recompile the executables. If you move the folders and wish to recompile, you should first delete the created "build" folder before recompiling.

Help
========================================

```
Usage: KnotAli[options] [input file]
```

Read input file from cmdline; predict minimum free energy and optimum structure using the time- and space-efficient MFE RNA folding algorithm for each sequence in the alignment.

```
  -h, --help             Print help and exit
  -V, --version          Print version and exit
  -v, --verbose          Turn on verbose
  -i, --input-type       Specify input file type (CLUSTAL or FASTA, base is FASTA)
  -o, --output-file      Specify output file
```

The input sequence is read from standard input, unless it is
given on the command line.

#### Example:
    assume you are in the source directory
    ./build/src/KnotAli myinputfile.afa
    ./build/src/KnotAli -o outputfile.afa myinputfile.afa
   
