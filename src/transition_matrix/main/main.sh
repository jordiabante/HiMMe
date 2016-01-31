#/usr/bin/env bash

rm -rf out

set -x
../transition_matrix.sh -k 1 -d out -- test.fa.gz
../transition_matrix.sh -k 2 -d out -- test.fa.gz
