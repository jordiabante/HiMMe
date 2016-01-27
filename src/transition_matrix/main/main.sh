#/usr/bin/env bash

rm -rf out

../transition_matrix.sh -t 4 -k 1 -d out -- test.fa
