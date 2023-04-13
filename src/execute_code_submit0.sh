#!/bin/bash
#PBS -l nodes=1:ppn=24
#PBS 
#PBS -j oe
export MKL_NUM_THREADS=24
cd $PBS_O_WORKDIR
cat $PBS_NODEFILE


./lri.out 256 1


#g++ sparsetest.cpp -o sparsetestclc -O3 -std=c++11 -larmadillo -I/home/opt/armadillo-9.600.5/usr/include/-L/home/opt/armadillo-9.600.5/usr/lib64/-Wl,-rpath=/home/opt/armadillo-9.600.5/usr/lib64/









