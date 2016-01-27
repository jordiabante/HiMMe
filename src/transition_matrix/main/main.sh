#/usr/bin/env bash

rm -rf out

set -x
../transition_matrix.sh -k 1 -d out -- test.fa
../transition_matrix.sh -k 2 -d out -- test.fa
