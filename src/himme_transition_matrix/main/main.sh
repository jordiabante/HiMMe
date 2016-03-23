#/usr/bin/env bash

rm -rf out

set -x
../himme_transition_matrix.sh -k 1 -d out -- test.fa.gz
../himme_transition_matrix.sh -k 2 -d out -- test.fa.gz
