#/usr/bin/env bash

rm -rf out

set -x
../himme_scoring.sh -t 1 -k 1 -d out \
    grch37_cdna_tm1.txt.gz emission_test.txt.gz input.fa
../himme_scoring.sh -t 1 -k 3 -d out \
    grch37_cdna_tm3.txt.gz emission_test.txt.gz input.fa
