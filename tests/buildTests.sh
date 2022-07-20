#!/bin/bash
#this short script builds darcyTester automaticaly

#this line sets gcc to default compiler on cluster, if you are not on the cluster comment this line out
module load gcc

rm -rf ./build
mkdir build
cd ./build
cmake ..

cmake --build .
