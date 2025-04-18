# Welcome to ABE

## 🚀 Welcome to the repo of ABE

The source code for the WISE paper under review: *Finding Top-k Keywords-Aware Optimal Routes: An Aggregation-Based Expansion Approach.*

## 🏠 Overview  

we study the k-KAOR problem and introduce an Aggregated-Based Expansion (ABE) algorithm that can support returning the top-k keywords-aware optimal routes.

## 1️⃣ Requirement

```C++
C++ 17
cmake 3.22
boost 1.8
metis 
```

## 1️⃣ Index File

Due to the large size of the index file (>1GB), we do not directly provide the index file, but instead provide the source code for building the index. This can also be used to build indexes locally.

BuildH2HIndex.cpp has built an H2H index for conducting shortest distance queries.
BuildIGTreeIndex.cpp has built an IGTree index.

Please note that these two files are separate C++ projects that can be run directly

##  3️⃣ How to run?

```cmake
mkdir build
cd build
cmake..
make
./ABE
```

