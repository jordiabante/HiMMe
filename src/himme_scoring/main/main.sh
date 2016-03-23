#/usr/bin/env bash

rm -rf out

set -x
../himme_scoring.sh -t 3 -k 1 -d out \
    transition_matrix.txt.gz emission_matrix.txt.gz input.fa
