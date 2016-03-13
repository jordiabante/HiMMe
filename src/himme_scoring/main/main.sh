#/usr/bin/env bash

rm -rf out

set -x
../himme_scoring.sh -k 3 -d out \
    grch37_cdna_tm3.txt.gz emission_test.txt.gz input.fa
